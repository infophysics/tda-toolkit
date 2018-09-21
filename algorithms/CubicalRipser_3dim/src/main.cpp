int main(int argc, char** argv){

	const char* filename = nullptr;
	string output_filename = "answer_3dim.diagram"; //default name
	file_format format = DIPHA;
	calculation_method method = LINKFIND;
	double threshold = 99999;
	bool print = false;

	for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--threshold") {
			string parameter = string(argv[++i]);
			size_t next_pos;
			threshold = stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--format") {
			string parameter = string(argv[++i]);
			if (parameter == "dipha") {
				format = DIPHA;
			} else if (parameter == "perseus") {
				format = PERSEUS;
			} else {
				print_usage_and_exit(-1);
			}
		} else if(arg == "--method") {
			string parameter = string(argv[++i]);
			if (parameter == "link_find") {
				method = LINKFIND;
			} else if (parameter == "compute_pairs") {
				method = COMPUTEPAIRS;
			} else {
				print_usage_and_exit(-1);
			}
		} else if (arg == "--output") {
			output_filename = string(argv[++i]);
		} else if (arg == "--print"){
			print = true;
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}
	
	ripser = CubicalRipser3D();
	ripser.ComputeBarcode(filename, output_filename, format, method, threshold, print);
	return 0;
}