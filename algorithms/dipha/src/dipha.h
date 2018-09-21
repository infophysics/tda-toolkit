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

class Dipha{

	public:
		Dipha();
		virtual ~Dipha();
		void print_help_and_exit();
		void parseCommandLine(int argc, char** argv, int type, bool& benchmark, bool& dualize, int64_t& upper_dim, std::string& input_filename, 
                        std::string& output_filename);
		void compute(const std::string& input_filename, bool dualize, int64_t upper_dim, const std::string& output_filename);
		void createPhatFiltration(const std::string& input_filename, bool dualize, int64_t upper_dim, const std::string& output_filename);
		void createSparseRepresentation(const std::string& input_filename, double upper_value, const std::string& output_filename);
};
