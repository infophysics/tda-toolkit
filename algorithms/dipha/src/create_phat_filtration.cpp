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

#include <dipha/includes.h>

int main(int argc, char** argv)
{
  // mandatory MPI initilization call
  MPI_Init(&argc, &argv);

  if (dipha::mpi_utils::get_num_processes() != 1)
  {
    dipha::mpi_utils::error_printer_if_root() << "Error: supports only one process." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  std::string input_filename; // name of file that contains the weighted cell complex
  std::string output_filename; // name of file that will contain the PHAT filtration
  bool dualize = false; // primal / dual computation toggle
  int64_t upper_dim = std::numeric_limits< int64_t >::max();
  dip = Dipha();
  dip.parseCommandLine(argc, argv, 1, dualize, upper_dim, input_filename, output_filename);

  switch (dipha::file_types::get_file_type(input_filename))
  {
  case dipha::file_types::IMAGE_DATA:
    if (upper_dim != std::numeric_limits< int64_t >::max())
    {
      dipha::mpi_utils::error_printer_if_root() << "upper_dim not supported for this input type IMAGE_DATA" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    dip.createPhatFiltration< dipha::inputs::weighted_cubical_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::WEIGHTED_BOUNDARY_MATRIX:
    if (upper_dim != std::numeric_limits< int64_t >::max())
    {
      dipha::mpi_utils::error_printer_if_root() << "upper_dim not supported for this input type WEIGHTED_BOUNDARY_MATRIX" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    dip.createPhatFiltration< dipha::inputs::weighted_explicit_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::DISTANCE_MATRIX:
    dip.createPhatFiltration< dipha::inputs::full_rips_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::SPARSE_DISTANCE_MATRIX:
    dip.createPhatFiltration< dipha::inputs::sparse_rips_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  default:
    dipha::mpi_utils::error_printer_if_root() << "Unknown complex type in DIPHA file" << input_filename << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  MPI_Finalize();
}
