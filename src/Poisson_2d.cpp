// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

#include <stdio.h>
#include <functional>
#include <list>
#include <chrono>
#include <ctime>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <tuple>
#include <iomanip>

#include <boost/random.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>


using namespace std;

#ifndef POISSON_2D_H
#define POISSON_2D_H

struct Cell_2d
{

  vector<double> coords_x;
  vector<double> coords_y;
  vector<double> death_rates;

  Cell_2d() {}
};

struct Grid_2d
{
  std::vector<Cell_2d> cells;
  std::vector<double> cell_death_rates;
  std::vector<int> cell_population;

  double area_length_x;
  double area_length_y;

  int cell_count_x;
  int cell_count_y;
  
  double b, d, dd;

  int seed; //?

  boost::random::lagged_fibonacci2281 rng;

  double initial_density;

  int total_population;
  double total_death_rate;

  double time;
  int event_count;

  std::vector<double> death_kernel_y;  //?
  std::vector<double> birth_kernel_y;  //?

  double spline_precision;  //?

  double death_cutoff_r;  //?
  double birth_cutoff_r;  //?

  double death_step;  //?
  double birth_step;  //?
  double birth_reverse_cdf_step;  //?

  int death_spline_nodes;  //?
  int birth_spline_nodes;  //?
  int birth_reverse_cdf_nodes;  //?

  boost::math::cubic_b_spline<double> death_kernel_spline;  //?
  boost::math::cubic_b_spline<double> birth_kernel_spline;  //?

  boost::math::cubic_b_spline<double> birth_reverse_cdf_spline;  //?

  int cull_x;
  int cull_y;

  std::tuple<std::tuple<double, double>, int> last_event;

  Cell_2d &cell_at(int i, int j)
  {
    return cells[cell_count_x*i + j];
  }

  double &cell_death_rate_at(int i, int j)
  {
    return cell_death_rates[cell_count_x*i + j];
  }

  int &cell_population_at(int i, int j)
  {
    return cell_population[cell_count_x*i + j];
  }

  vector<double> get_coords_x_at_cell(int i, int j)
  {
    return cells[cell_count_x*i + j].coords_x;
  }
  
  vector<double> get_coords_y_at_cell(int i, int j)
  {
    return cells[cell_count_x*i + j].coords_y;
  }

  vector<double> get_death_rates_at_cell(int i, int j)
  {
    return cells[cell_count_x*i + j].death_rates;
  }

  vector<double> get_all_coords_x()
  {
    vector<double> result;
    for (auto cell : cells)
    {
      if (cell.coords_x.size() != 0)
        result.insert(result.end(), cell.coords_x.begin(), cell.coords_x.end());
    }
    return result;
  }

  vector<double> get_all_coords_y()
  {
    vector<double> result;
    for (auto cell : cells)
    {
      if (cell.coords_x.size() != 0)
        result.insert(result.end(), cell.coords_y.begin(), cell.coords_y.end());
    }
    return result;
  }
  
  vector<double> get_all_death_rates()
  {
    vector<double> result;
    for (auto cell : cells)
    {
      if (cell.coords_x.size() != 0)
        result.insert(result.end(), cell.death_rates.begin(), cell.death_rates.end());
    }
    return result;
  }

  void Initialize_death_rates()
  {

    for (int i = 0; i < cell_count_x * cell_count_y; i++)
    {
      cells.push_back(Cell_2d());
      cell_death_rates.push_back(0);
      cell_population.push_back(0);
    }

    //Spawn all speciments
    total_population = static_cast<int>(ceil(area_length_x*area_length_y*initial_density)); //initial population at t=0

    {
      double x_coord, y_coord;
      int i, j;
      for (int k = 0; k < total_population; k++)
      {
        x_coord = boost::uniform_real<>(0, area_length_x)(rng);
        y_coord = boost::uniform_real<>(0, area_length_y)(rng);
        
        i = static_cast<int>(floor(x_coord * cell_count_x / area_length_x));
        j = static_cast<int>(floor(y_coord * cell_count_y / area_length_y));
        
        if (i == cell_count_x)
          i--;
        if (j == cell_count_y)
          j--;

        cell_at(i, j).coords_x.push_back(x_coord);
        cell_at(i, j).coords_y.push_back(y_coord);
        
        cell_at(i, j).death_rates.push_back(d);
        cell_death_rate_at(i, j) += d;
        total_death_rate += d;

        cell_population_at(i, j)++;
      }
    }

    for (int i = 0; i < cell_count_x; i++)
    {
      for (int j = 0; j < cell_count_y; j++)
      {
        for (int k = 0; k < cell_population_at(i, j); k++)
        {
          for (int n = max(0, i - cull_x); n < min(cell_count_x, i + cull_x + 1); n++)
          {
            for (int m = max(0, j - cull_y); m < min(cell_count_y, j + cull_y + 1); m++)
            {
              for (int p = 0; p < cell_population_at(n, m); p++)
              {
                if (i == n && j == m && k == p)
                  continue; // same speciment
    
                //Distance between k-th speciment in (i, j) cell and p-th speciment in (n, m) cell
    
                double distance_x = abs(cell_at(i, j).coords_x[k] - cell_at(n, m).coords_x[p]);
                double distance_y = abs(cell_at(i, j).coords_y[k] - cell_at(n, m).coords_y[p]);
                
                double distance = sqrt(distance_x*distance_x + distance_y*distance_y);

                if (distance > death_cutoff_r)
                  continue; //Too far to interact
    
                double interaction = dd * death_kernel_spline(distance);
    
                cell_at(i, j).death_rates[k] += interaction;
                cell_death_rate_at(i, j) += interaction;
                total_death_rate += interaction;
              }
            }
          }
        }
      }
    }
  }

  void kill_random()
  {

    if (total_population == 1)
    {
      total_population--;
      return;
    }

    int cell_death_index = boost::random::discrete_distribution<>(cell_death_rates)(rng);
    int in_cell_death_index = boost::random::discrete_distribution<>(cells[cell_death_index].death_rates)(rng);

    Cell_2d &death_cell = cells[cell_death_index];

    last_event = make_tuple<>(make_tuple<>(death_cell.coords_x[in_cell_death_index], death_cell.coords_y[in_cell_death_index]), -1);

    int cell_death_x = cell_death_index / cell_count_x;
    int cell_death_y = cell_death_index % cell_count_x;
    
    for (int i = max(0, cell_death_x - cull_x); i < min(cell_count_x, cell_death_x + cull_x + 1); i++)
    {
      for (int j = max(0, cell_death_y - cull_y); j < min(cell_count_y, cell_death_y + cull_y + 1); j++)
      {
        for (int k = 0; k < cell_population_at(i, j); k++)
        {
          if (i == cell_death_x && j == cell_death_y && k == in_cell_death_index)
            continue;

          double distance_x = abs(cell_at(i, j).coords_x[k] - death_cell.coords_x[in_cell_death_index]);
          double distance_y = abs(cell_at(i, j).coords_y[k] - death_cell.coords_y[in_cell_death_index]);
          
          double distance = distance_x*distance_x + distance_y*distance_y;
  
          if (distance > death_cutoff_r*death_cutoff_r)
            continue; //Too far to interact
  
          double interaction = dd * death_kernel_spline(distance); //?
  
          cell_at(i, j).death_rates[k] -= interaction;
          //ignore dying speciment death rates since it is to be deleted
  
          cell_death_rate_at(i, j) -= interaction;
          cell_death_rate_at(cell_death_x, cell_death_y) -= interaction;
  
          total_death_rate -= 2 * interaction;
        }
      }
    }
    //remove dead speciment
    cell_death_rates[cell_death_index] -= d;
    total_death_rate -= d;

    if (abs(cell_death_rates[cell_death_index]) < 1e-10)
    {
      cell_death_rates[cell_death_index] = 0;
    }

    cell_population[cell_death_index]--;
    total_population--;

    //swap dead and last
    death_cell.death_rates[in_cell_death_index] = death_cell.death_rates[death_cell.death_rates.size() - 1];
    death_cell.coords_x[in_cell_death_index] = death_cell.coords_x[death_cell.coords_x.size() - 1];
    death_cell.coords_y[in_cell_death_index] = death_cell.coords_y[death_cell.coords_y.size() - 1];
    
    death_cell.death_rates.erase(death_cell.death_rates.end() - 1);
    death_cell.coords_x.erase(death_cell.coords_x.end() - 1);
    death_cell.coords_y.erase(death_cell.coords_y.end() - 1);
  }

  void spawn_random()
  {
    int cell_index = boost::random::discrete_distribution<>(cell_population)(rng);

    int event_index = boost::random::uniform_smallint<>(0, cell_population[cell_index] - 1)(rng);

    Cell_2d &parent_cell = cells[cell_index];
    
    double Theta = boost::uniform_real<>(0, 2*M_PI)(rng);
    double x_coord_new = parent_cell.coords_x[event_index] + birth_reverse_cdf_spline(boost::uniform_real<>(0, 1)(rng))*cos(Theta);
    double y_coord_new = parent_cell.coords_y[event_index] + birth_reverse_cdf_spline(boost::uniform_real<>(0, 1)(rng))*sin(Theta);

    
    if (x_coord_new < 0 || x_coord_new > area_length_x || y_coord_new < 0 || y_coord_new > area_length_y)
    {
      last_event = make_tuple<>(make_tuple<>(x_coord_new, y_coord_new), 0);
      //Speciment failed to spawn and died outside area boundaries
    }
    else
    {
      last_event = make_tuple<>(make_tuple<>(x_coord_new, y_coord_new), 1);

      int new_i = static_cast<int>(floor(x_coord_new * cell_count_x / area_length_x));
      int new_j = static_cast<int>(floor(y_coord_new * cell_count_y / area_length_y));
      
      if (new_i == cell_count_x)
        new_i--;

      if (new_j == cell_count_y)
        new_j--;
      
      //New speciment is added to the end of vector

      cell_at(new_i, new_j).coords_x.push_back(x_coord_new);
      cell_at(new_i, new_j).coords_y.push_back(y_coord_new);
      cell_at(new_i, new_j).death_rates.push_back(d);

      cell_death_rate_at(new_i, new_j) += d;
      total_death_rate += d;

      cell_population_at(new_i, new_j)++;
      total_population++;

      for (int i = max(0, new_i - cull_x); i < min(cell_count_x, new_i + cull_x + 1); i++)
      {
        for (int j = max(0, new_j - cull_y); j < min(cell_count_y, new_j + cull_y + 1); j++)
        {
          for (int k = 0; k < cell_population_at(i, j); k++)
          {
            if (i == new_i && j == new_j && k == cell_population_at(new_i, new_j) - 1)
              continue;
  
            double distance_x = abs(cell_at(i, j).coords_x[k] - x_coord_new);
            double distance_y = abs(cell_at(i, j).coords_y[k] - y_coord_new);
            
            double distance = distance_x*distance_x + distance_y*distance_y;
  
            if (distance > death_cutoff_r*death_cutoff_r)
              continue; //Too far to interact
  
            double interaction = dd * death_kernel_spline(distance); //?
  
            cell_at(i, j).death_rates[k] += interaction;
            cell_at(new_i, new_j).death_rates[cell_population_at(new_i, new_j) - 1] += interaction;
  
            cell_death_rate_at(i, j) += interaction;
            cell_death_rate_at(new_i, new_j) += interaction;
  
            total_death_rate += 2 * interaction;
          }
        }
      }
    }
  }

  void make_event()
  {
    if (total_population == 0)
      return;
    event_count++;
    time += boost::random::exponential_distribution<>(total_population * b + total_death_rate)(rng);
    //Rolling event according to global birth \ death rate
    if (boost::random::bernoulli_distribution<>(total_population * b / (total_population * b + total_death_rate))(rng) == 0)
    {
      kill_random();
    }
    else
    {
      spawn_random();
    }
  }
  void run_events(int events)
  {
    if (events > 0)
    {
      for (int i = 0; i < events; i++)
      {
        make_event();
      }
    }
  }

  void run_for(double time)
  {
    if (time > 0.0)
    {
      double time0 = this->time;
      while (this->time < time0 + time)
      {
        make_event();
        if (total_population == 0)
          return;
      }
    }
  }

  double get_birth_spline_value(double at) //?
  {
    return birth_kernel_spline(at);
  }

  double get_death_spline_value(double at) //?
  {
    return death_kernel_spline(at);
  }

  double get_birth_reverse_cdf_spline_value(double at) //?
  {
    return birth_reverse_cdf_spline(at);
  }

  
  Grid_2d(Rcpp::List params) : time(), event_count(), cells(), cell_death_rates(), cell_population(),
                               death_kernel_spline(), birth_kernel_spline(), birth_reverse_cdf_spline()
  {

    //Parse parameters

    Rcpp::Environment base("package:base");

    // Make function callable from C++
    Rcpp::Function print = base["print"];

    area_length_x = Rcpp::as<double>(params["area_length_x"]);
    area_length_y = Rcpp::as<double>(params["area_length_y"]);
    cell_count_x = Rcpp::as<int>(params["cell_count_x"]);
    cell_count_y = Rcpp::as<int>(params["cell_count_y"]);

    b = Rcpp::as<double>(params["b"]);
    d = Rcpp::as<double>(params["d"]);
    dd = Rcpp::as<double>(params["dd"]);

    seed = Rcpp::as<int>(params["seed"]);
    rng = boost::random::lagged_fibonacci2281(uint32_t(seed));

    initial_density = Rcpp::as<double>(params["init_density"]);

    //?
    
    death_kernel_y = Rcpp::as<vector<double>>(params["death_kernel_y"]);
    death_cutoff_r = Rcpp::as<double>(params["death_kernel_r"]);
    death_spline_nodes = death_kernel_y.size();
    death_step = death_cutoff_r/(death_spline_nodes-1);

    birth_kernel_y = Rcpp::as<vector<double>>(params["birth_kernel_y"]);
    birth_cutoff_r = Rcpp::as<double>(params["birth_kernel_r"]);
    birth_spline_nodes = birth_kernel_y.size();
    birth_step = birth_cutoff_r/(birth_spline_nodes-1);

    spline_precision = Rcpp::as<double>(params["spline_precision"]);

    using boost::math::cubic_b_spline;
    //Build death spline, ensure 0 derivative at 0 (symmetric) and endpoint (expected no death interaction further)
    death_kernel_spline=cubic_b_spline<double>(death_kernel_y.begin(), death_kernel_y.end(), 0, death_step, 0, 0);
    //trim_spline(death_kernel_x, death_kernel_y, death_kernel_spline, death_spline_nodes, death_cutoff_r, spline_precision);

    //Calculate amount of cells to check around for death interaction

    cull_x = max(static_cast<int>(ceil(death_cutoff_r / (area_length_x / cell_count_x))), 3);
    cull_y = max(static_cast<int>(ceil(death_cutoff_r / (area_length_y / cell_count_y))), 3);
    
    //Build birth spline
    birth_kernel_spline=cubic_b_spline<double>(birth_kernel_y.begin(), birth_kernel_y.end(), 0, birth_step, 0, 0);
    //trim_spline(birth_kernel_x, birth_kernel_y, birth_kernel_spline, birth_spline_nodes, birth_cutoff_r, spline_precision);
    //Calculate reverse CDF for birth spline

    vector<double> x_quantile_1d_array(birth_spline_nodes);
    vector<double> y_quantile_1d_array(birth_spline_nodes);

    using boost::math::quadrature::trapezoidal;
    double approx_const = trapezoidal([=](double y) {return birth_kernel_spline(y);}, 0.0, birth_cutoff_r);

    using boost::math::tools::newton_raphson_iterate;

    for (int i = 0; i < birth_spline_nodes; i++)
    {
      x_quantile_1d_array[i] = (double)i / (birth_spline_nodes - 1);
      y_quantile_1d_array[i] =
          newton_raphson_iterate(
              [=](double y) { return make_tuple(
                                  trapezoidal([=](double z) {return birth_kernel_spline(z);}, 0.0, y) / approx_const - x_quantile_1d_array[i],
                                  birth_kernel_spline(y) / approx_const); },
              1e-10, 0.0, birth_cutoff_r, numeric_limits<double>::digits);
    }
    int i = 0;
    while (y_quantile_1d_array[i] < 1e-300)
    {
      i++;
    }

    vector<double> y_quantile_1d_array_temp(birth_spline_nodes - i);

    for (int j = 0; j < birth_spline_nodes - i; j++)
    {
      y_quantile_1d_array_temp[j] = y_quantile_1d_array[i + j];
    }
    
    birth_reverse_cdf_step = 1.0/(y_quantile_1d_array_temp.size()-1);
    birth_reverse_cdf_nodes = y_quantile_1d_array_temp.size();
    // Extrapolate last quantile element
    double right_derivative = (y_quantile_1d_array_temp[y_quantile_1d_array_temp.size()-2] - y_quantile_1d_array_temp[y_quantile_1d_array_temp.size()-3])/birth_reverse_cdf_step;
    
    y_quantile_1d_array_temp[y_quantile_1d_array_temp.size()-1] = y_quantile_1d_array_temp[y_quantile_1d_array_temp.size()-2] + right_derivative * birth_reverse_cdf_step;
      
    //Ensure correct derivative at 0, equal to 1/2*birth_kernel(0)
    birth_reverse_cdf_spline = cubic_b_spline<double>(y_quantile_1d_array_temp.begin(), y_quantile_1d_array_temp.end(), 0, birth_reverse_cdf_step, 
                                                      0.5/birth_kernel_spline(0), 2*right_derivative);

    //?
    
    //Spawn speciments and calculate death rates
    Initialize_death_rates();
    total_death_rate = accumulate(cell_death_rates.begin(), cell_death_rates.end(), 0.0);
  }
};

RCPP_MODULE(poisson_2d_module)
{
  using namespace Rcpp;

  class_<Grid_2d>("poisson_2d")
      .constructor<List>()
      .field_readonly("area_length_x", &Grid_2d::area_length_x)
      .field_readonly("area_length_y", &Grid_2d::area_length_y)
      .field_readonly("cell_count_x", &Grid_2d::cell_count_x)
      .field_readonly("cell_count_y", &Grid_2d::cell_count_y)
  
      .field_readonly("b", &Grid_2d::b)
      .field_readonly("d", &Grid_2d::d)
      .field_readonly("dd", &Grid_2d::dd)

      .field_readonly("seed", &Grid_2d::seed)
      .field_readonly("initial_density", &Grid_2d::initial_density)

      .field_readonly("death_kernel_y", &Grid_2d::death_kernel_y)
      .field_readonly("death_cutoff_r", &Grid_2d::death_cutoff_r)
      .field_readonly("death_spline_nodes", &Grid_2d::death_spline_nodes)
      .field_readonly("death_step", &Grid_2d::death_step)

      .field_readonly("birth_kernel_y", &Grid_2d::birth_kernel_y)
      .field_readonly("birth_cutoff_r", &Grid_2d::birth_cutoff_r)
      .field_readonly("birth_spline_nodes", &Grid_2d::birth_spline_nodes)
      .field_readonly("birth_step", &Grid_2d::birth_step)

      .field_readonly("birth_reverse_cdf_nodes", &Grid_2d::birth_reverse_cdf_nodes)
      .field_readonly("birth_reverse_cdf_step", &Grid_2d::birth_reverse_cdf_step)

      .field_readonly("spline_precision", &Grid_2d::spline_precision)

      .field_readonly("cell_death_rates", &Grid_2d::cell_death_rates)
      .field_readonly("cell_population", &Grid_2d::cell_population)

      .method("get_all_x_coordinates", &Grid_2d::get_all_coords_x)
      .method("get_all_y_coordinates", &Grid_2d::get_all_coords_y)
      .method("get_all_death_rates", &Grid_2d::get_all_death_rates)

      .method("get_x_coordinates_in_cell", &Grid_2d::get_coords_x_at_cell)
      .method("get_y_coordinates_in_cell", &Grid_2d::get_coords_y_at_cell)
      .method("get_death_rates_in_cell", &Grid_2d::get_death_rates_at_cell)

      .method("birth_spline_at", &Grid_2d::get_birth_spline_value)
      .method("death_spline_at", &Grid_2d::get_death_spline_value)
      .method("birth_reverse_cdf_spline_at", &Grid_2d::get_birth_reverse_cdf_spline_value)

      .method("make_event", &Grid_2d::make_event)
      .method("run_events", &Grid_2d::run_events)
      .method("run_for", &Grid_2d::run_for)

      .field_readonly("total_population", &Grid_2d::total_population)
      .field_readonly("total_death_rate", &Grid_2d::total_death_rate)
      .field_readonly("events", &Grid_2d::event_count)
      .field_readonly("time", &Grid_2d::time);
}

#endif
