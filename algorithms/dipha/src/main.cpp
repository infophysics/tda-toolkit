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

  // take time at beggining of execution
  double time_at_start = MPI_Wtime();

  std::string input_filename; // name of file that contains the weighted cell complex
  std::string output_filename; // name of file that will contain the persistence diagram
  bool benchmark = false; // print timings / info
  bool dualize = false; // primal / dual computation toggle
  int64_t upper_dim = std::numeric_limits< int64_t >::max();
  dip = Dipha();
  dip.parseCommandLine(argc, argv, 0, benchmark, dualize, upper_dim, input_filename, output_filename);

  if (benchmark)
  {
    dipha::globals::benchmark = true;

    dipha::mpi_utils::cout_if_root() << std::endl << "Input filename: " << std::endl << input_filename << std::endl;

    dipha::mpi_utils::cout_if_root() << std::endl << "upper_dim: " << upper_dim << std::endl;

    dipha::mpi_utils::cout_if_root() << std::endl << "Number of processes used: " << std::endl << dipha::mpi_utils::get_num_processes() 
                                     << std::endl;

    dipha::mpi_utils::cout_if_root() << std::endl << "Detailed information for rank 0:" << std::endl;

    dipha::mpi_utils::cout_if_root() << std::setw(11) << "time" << std::setw(13) << "prior mem" << std::setw(13) << "peak mem" 
                                     << std::setw(13) << "bytes recv" << std::endl;
  }

  switch (dipha::file_types::get_file_type(input_filename))
  {
  case dipha::file_types::IMAGE_DATA:
    dip.compute< dipha::inputs::weighted_cubical_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::WEIGHTED_BOUNDARY_MATRIX:
    dip.compute< dipha::inputs::weighted_explicit_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::DISTANCE_MATRIX:
    dip.compute< dipha::inputs::full_rips_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::SPARSE_DISTANCE_MATRIX:
    dip.compute< dipha::inputs::sparse_rips_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  default:
    dipha::mpi_utils::error_printer_if_root() << "Unknown complex type in DIPHA file" << input_filename << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  if (benchmark)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    dipha::mpi_utils::cout_if_root() << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) << std::setprecision(1);
    dipha::mpi_utils::cout_if_root() << std::endl << "Overall running time in seconds: " << std::endl;
    dipha::mpi_utils::cout_if_root() << std::setprecision(1) << MPI_Wtime() - time_at_start << std::endl;

    dipha::mpi_utils::cout_if_root() << std::endl << "Reduction kernel running time in seconds: " << std::endl;
    dipha::mpi_utils::cout_if_root() << std::setprecision(1) << dipha::globals::reduction_kernel_running_time << std::endl;

    int64_t peak_mem = getPeakRSS() >> 20;
    std::vector< int64_t > peak_mem_per_rank(dipha::mpi_utils::get_num_processes());
    MPI_Gather(&peak_mem, 1, MPI_LONG_LONG, peak_mem_per_rank.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    dipha::mpi_utils::cout_if_root() << std::endl << "Overall peak mem in GB of all ranks: " << std::endl;
    dipha::mpi_utils::cout_if_root() << (double)*std::max_element(peak_mem_per_rank.begin(), peak_mem_per_rank.end()) / 1024.0 
                                     << std::endl;

    dipha::mpi_utils::cout_if_root() << std::endl << "Individual peak mem in GB of per rank: " << std::endl;
    for (int64_t peak_mem : peak_mem_per_rank)
    {
      dipha::mpi_utils::cout_if_root() << (double)peak_mem / 1024.0 << std::endl;
    }

    std::vector< int64_t > bytes_received_per_rank(dipha::mpi_utils::get_num_processes());
    MPI_Gather(&dipha::globals::bytes_received, 1, MPI_LONG_LONG, bytes_received_per_rank.data(), 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    dipha::mpi_utils::cout_if_root() << std::endl << "Maximal communication traffic (without sorting) in GB between any pair of nodes:" 
                                     << std::endl;
    const int64_t max_num_MB_traffic = *std::max_element(bytes_received_per_rank.begin(), bytes_received_per_rank.end()) >> 20;
    dipha::mpi_utils::cout_if_root() << std::setprecision(1) << (double)(max_num_MB_traffic) / 1024.0 << std::endl;

    dipha::mpi_utils::cout_if_root() << std::endl << "Total communication traffic (without sorting) in GB between all pairs of nodes:" 
                                     << std::endl;
    const int64_t total_num_MB_traffic = std::accumulate(bytes_received_per_rank.begin(), bytes_received_per_rank.end(), 0LL) >> 20;
    dipha::mpi_utils::cout_if_root() << std::setprecision(1) << (double)(total_num_MB_traffic) / 1024.0 << std::endl;
  }

  MPI_Finalize();
}
