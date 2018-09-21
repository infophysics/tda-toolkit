/*  Copyright 2014 IST Austria

Contributed by: Jan Reininghaus

This file is part of DIPHA.

DIPHA is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DIPHA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with DIPHA.  If not, see <http://www.gnu.org/licenses/>. */

#include "../include/dipha/includes.h"

void Dipha::print_help_and_exit()
{
  std::cerr << "Usage: " << "dipha " << "[options] input_filename output_filename" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << std::endl;
  std::cerr << "--help    --  prints this screen" << std::endl;
  std::cerr << "--dual    --  use dualization" << std::endl;
  std::cerr << "--upper_dim N   --  maximal dimension to compute" << std::endl;
  std::cerr << "--benchmark --  prints timing info" << std::endl;
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void Dipha::parseCommandLine(int argc, char** argv, int type, bool& benchmark, bool& dualize, int64_t& upper_dim, std::string& input_filename, 
                        std::string& output_filename)
{
  if (type == 0){

	  if (argc < 3)
	    print_help_and_exit();

	  input_filename = argv[argc - 2];
	  output_filename = argv[argc - 1];

	  for (int idx = 1; idx < argc - 2; idx++)
	  {
	    const std::string option = argv[idx];
	    if (option == "--benchmark")
	    {
	      benchmark = true;
	    }
	    else if (option == "--help")
	    {
	      print_help_and_exit();
	    }
	    else if (option == "--dual")
	    {
	      dualize = true;
	    }
	    else if (option == "--upper_dim")
	    {
	      idx++;
	      if (idx >= argc - 2)
	        print_help_and_exit();
	      std::string parameter = std::string(argv[idx]);
	      size_t pos_last_digit;
	      upper_dim = std::stoll(parameter, &pos_last_digit);
	      if (pos_last_digit != parameter.size())
	        print_help_and_exit();
	    }
	    else print_help_and_exit();
	  }
	}
	else if (type == 1){
	if (argc < 3)
    print_help_and_exit();

  input_filename = argv[argc - 2];
  output_filename = argv[argc - 1];

  for (int idx = 1; idx < argc - 2; idx++)
  {
    const std::string option = argv[idx];
    if (option == "--help")
    {
      print_help_and_exit();
    }
    else if (option == "--dual")
    {
      dualize = true;
    }
    else if (option == "--upper_dim")
    {
      idx++;
      if (idx >= argc - 2)
        print_help_and_exit();
      std::string parameter = std::string(argv[idx]);
      size_t pos_last_digit;
      upper_dim = std::stoll(parameter, &pos_last_digit);
      if (pos_last_digit != parameter.size())
        print_help_and_exit();
    }
    else print_help_and_exit();
  }
	}
	else{
		if (argc < 3)
    print_help_and_exit();

  input_filename = argv[2];
  output_filename = argv[3];

  std::string parameter = std::string(argv[1]);
  size_t pos_last_digit;
  upper_value = std::stod(parameter, &pos_last_digit);
  if (pos_last_digit != parameter.size())
    print_help_and_exit();
	}

}

template< typename Complex >
void Dipha::compute(const std::string& input_filename, bool dualize, int64_t upper_dim, const std::string& output_filename)
{
  Complex complex;
  DIPHA_MACROS_BENCHMARK(complex.load_binary(input_filename, upper_dim); );
  if (dipha::globals::benchmark)
    dipha::mpi_utils::cout_if_root() << std::endl << "Number of cells in input: " << std::endl << complex.get_num_cells() << std::endl;
  dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
  dipha::data_structures::write_once_column_array reduced_columns;
  dipha::algorithms::compute_reduced_columns(complex, dualize, upper_dim, filtration_to_cell_map, reduced_columns);
  DIPHA_MACROS_BENCHMARK(dipha::outputs::save_persistence_diagram(output_filename, complex, filtration_to_cell_map, reduced_columns, 
                                                                  dualize, upper_dim); );
}
void Dipha::createPhatFiltration(const std::string& input_filename, bool dualize, int64_t upper_dim, const std::string& output_filename)
{
  Complex complex;
  complex.load_binary(input_filename, upper_dim);
  dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
  dipha::algorithms::get_filtration_to_cell_map(complex, dualize, filtration_to_cell_map);
  dipha::data_structures::distributed_vector< int64_t > cell_to_filtration_map;
  dipha::algorithms::get_cell_to_filtration_map(complex.get_num_cells(), filtration_to_cell_map, cell_to_filtration_map);

  const int64_t nr_columns = complex.get_num_cells();
  std::vector< std::vector< int64_t > > matrix(nr_columns);
  std::vector< int64_t > dims(nr_columns, -1);

  for (int64_t cur_dim = 0; cur_dim <= complex.get_max_dim(); cur_dim++)
  {
    dipha::data_structures::flat_column_stack unreduced_columns;
    dipha::algorithms::generate_unreduced_columns(complex, filtration_to_cell_map, cell_to_filtration_map, cur_dim, dualize, 
                                                  unreduced_columns);
    dipha::data_structures::heap_column col;
    while (!unreduced_columns.empty())
    {
      int64_t index;
      unreduced_columns.pop(index, col);
      std::sort(col.begin(), col.end());
      matrix[index] = col;
      dims[index] = dualize ? complex.get_max_dim() - cur_dim : cur_dim;
    }
  }

  std::ofstream output_stream(output_filename.c_str(), std::ios_base::binary | std::ios_base::out);
  output_stream.write((char*)&nr_columns, sizeof(int64_t));
  for (int64_t cur_col = 0; cur_col < nr_columns; cur_col++)
  {
    int64_t cur_dim = dims[cur_col];
    output_stream.write((char*)&cur_dim, sizeof(int64_t));
    int64_t cur_nr_rows = matrix[cur_col].size();
    output_stream.write((char*)&cur_nr_rows, sizeof(int64_t));
    for (int64_t cur_row_idx = 0; cur_row_idx < cur_nr_rows; cur_row_idx++)
    {
      int64_t cur_row = matrix[cur_col][cur_row_idx];
      output_stream.write((char*)&cur_row, sizeof(int64_t));
    }
  }

  output_stream.close();
}
void Dipha::createSparseRepresentation(const std::string& input_filename, double upper_value, const std::string& output_filename)
{
  dipha::inputs::full_rips_complex complex;
  complex.load_binary(input_filename, 1);

  std::vector< std::vector< std::pair<int64_t, double> > > sparse_matrix;

  int64_t matrix_size = complex.number_of_points();

  sparse_matrix.resize(matrix_size);

  for (int64_t rows = 0; rows < matrix_size; rows++)
  {

    for (int64_t cols = 0; cols < matrix_size; cols++)
    {

      double distance = complex.get_distance(rows, cols);

      if (rows != cols && distance <= 2 * upper_value)
      {
        sparse_matrix[rows].push_back(std::make_pair(cols, distance));
      }
    }
  }

  std::ofstream output_stream(output_filename.c_str(), std::ios_base::binary | std::ios_base::out);
  int64_t magic_number = 8067171840;
  output_stream.write((char*)&magic_number, sizeof(int64_t));
  int64_t file_format_id = dipha::file_types::SPARSE_DISTANCE_MATRIX;
  output_stream.write((char*)&file_format_id, sizeof(int64_t));


  output_stream.write((char*)&matrix_size, sizeof(int64_t));
  for (int64_t cur_row = 0; cur_row < matrix_size; cur_row++)
  {
    int64_t row_size = sparse_matrix[cur_row].size();
    output_stream.write((char*)&row_size, sizeof(int64_t));
  }
  for (int64_t cur_row = 0; cur_row < matrix_size; cur_row++)
  {
    int64_t row_size = sparse_matrix[cur_row].size();
    for (int64_t cur_col_idx = 0; cur_col_idx < row_size; cur_col_idx++)
    {
      output_stream.write((char*)&(sparse_matrix[cur_row][cur_col_idx].first), sizeof(int64_t));
    }
  }
  for (int64_t cur_row = 0; cur_row < matrix_size; cur_row++)
  {
    int64_t row_size = sparse_matrix[cur_row].size();
    for (int64_t cur_col_idx = 0; cur_col_idx < row_size; cur_col_idx++)
    {
      output_stream.write((char*)&(sparse_matrix[cur_row][cur_col_idx].second), sizeof(double));
    }
  }

  output_stream.close();
}

