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

  std::string input_filename; // name of file that contains the weighted cell complex
  std::string output_filename; // name of file that will contain the PHAT filtration
  double upper_value = std::numeric_limits< double >::max();
  dip = Dipha();
  dip.parseCommandLine(argc, argv, 2, upper_value, input_filename, output_filename);

  if (dipha::file_types::get_file_type(input_filename) != dipha::file_types::DISTANCE_MATRIX)
  {

    std::cerr << "Requires DIPHA Distance matrix file format for " << input_filename << std::endl;
    std::exit(-1);
  }
  dip.createSparseRepresentation(input_filename, upper_value, output_filename);

  MPI_Finalize();
}
