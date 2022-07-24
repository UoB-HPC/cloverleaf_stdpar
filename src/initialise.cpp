/*
 Crown Copyright 2012 AWE.

 This file is part of CloverLeaf.

 CloverLeaf is free software: you can redistribute it and/or modify it under 
 the terms of the GNU General Public License as published by the 
 Free Software Foundation, either version 3 of the License, or (at your option) 
 any later version.

 CloverLeaf is distributed in the hope that it will be useful, but 
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
 details.

 You should have received a copy of the GNU General Public License along with
 CloverLeaf. If not, see http://www.gnu.org/licenses/.
 */


//  @brief Top level initialisation routine
//  @author Wayne Gaudin
//  @details Checks for the user input and either invokes the input reader or
//  switches to the internal test problem. It processes the input and strips
//  comments before writing a final input file.
//  It then calls the start routine.

#include "initialise.h"
#include "report.h"
#include "version.h"
#include "read_input.h"
#include "start.h"

#include <fstream>
#include <iomanip>

extern std::ostream g_out;
std::ofstream of;


struct RunConfig {
	std::string file;
	size_t deviceIdx;
};


void printSimple(size_t num_devices) {
	for (size_t j = 0; j < num_devices; ++j) {
		std::cout << std::setw(3) << j << ". offload device #" << j << std::endl;
	}
}

void printHelp(const std::string &name) {
	std::cout << std::endl;
	std::cout << "Usage: " << name << " [OPTIONS]\n\n"
	          << "Options:\n"
	          << "  -h  --help               Print the message\n"
	          << "      --list               List available devices\n"
	          << "      --no-target          Use OMP fallback\n"
	          << "      --device <INDEX>     Select device at INDEX from output of --list\n"
	          << "      --file               Custom clover.in file (defaults to clover.in if unspecified)\n"
	          << std::endl;
}

RunConfig parseArgs(const size_t num_devices,
                    const std::vector<std::string> &args) {

	const auto readParam = [&args](size_t current, const std::string &emptyMessage, auto map) {
		if (current + 1 < args.size()) {
			return map(args[current + 1]);
		} else {
			std::cerr << emptyMessage << std::endl;
			printHelp(args[0]);
			std::exit(EXIT_FAILURE);
		}
	};

	auto config = RunConfig{"clover.in", 0};
	for (size_t i = 0; i < args.size(); ++i) {
		const auto arg = args[i];

		if (arg == "--help" || arg == "-h") {
			printHelp(args[0]);
			std::exit(EXIT_SUCCESS);
		} else if (arg == "--list") {
			printSimple(num_devices);
			std::exit(EXIT_SUCCESS);
		} else if (arg == "--no-target") {
			config.deviceIdx = -1;
		} else if (arg == "--device") {
			readParam(i, "--device specified but no index was given", [&](const auto &param) {
				auto selected = std::stoul(param);
				if (selected < 0 || selected >= num_devices) {
					std::cerr << "bad device index `" << param << "`" << std::endl;
					std::exit(EXIT_FAILURE);
				}
				config.deviceIdx = selected;
			});
		} else if (arg == "--file") {
			readParam(i, "--file specified but no file was given", [&config](const auto &param) {
				config.file = param;
			});
		}
	}
	return config;
}

std::unique_ptr<global_variables>
initialise(parallel_ &parallel, const std::vector<std::string> &args) {

	global_config config;

	if (parallel.boss) {
		of.open("clover.out");
		if (!of.is_open())
			report_error((char *) "initialise", (char *) "Error opening clover.out file.");

		g_out.rdbuf(of.rdbuf());

	} else {
		g_out.rdbuf(std::cout.rdbuf());
	}

	if (parallel.boss) {
		g_out << "Clover Version " << g_version << std::endl
		      << "Kokkos Version" << std::endl
		      << "Task Count " << parallel.max_task << std::endl
		      << std::endl;

		std::cout << "Output file clover.out opened. All output will go there." << std::endl;
	}

	clover_barrier();

    // XXX stdpar doesn't support listing devices
	auto num_devices = 1;
	if (num_devices == 0) {
		std::cout << "No offload devices available" << std::endl;
	} else {
		std::cout << "Detected offload devices:" << std::endl;
		printSimple(num_devices);
	}

	auto runConfig = parseArgs(num_devices, args);
	auto file = runConfig.file;
	auto selectedDevice = runConfig.deviceIdx;
	std::cout << "Using OMP device: " << selectedDevice << std::endl;

	std::ifstream g_in;
	if (parallel.boss) {
		g_out << "Clover will run from the following input:-" << std::endl
		      << std::endl;

		if (!args.empty()) {
			std::cout << "Args:";
			for (const auto &arg : args) std::cout << " " << arg;
			std::cout << std::endl;
		}

		g_in.open(file);
		std::cout << "Using input: `" << file << "`" << std::endl;

		if (!g_in.good()) {
			std::cerr << "Unable to open file: `" << file << "`, using defaults" << std::endl;
			g_in.close();
			std::ofstream out_unit("clover.in");
			out_unit
					<< "*clover" << std::endl
					<< " state 1 density=0.2 energy=1.0" << std::endl
					<< " state 2 density=1.0 energy=2.5 geometry=rectangle xmin=0.0 xmax=5.0 ymin=0.0 ymax=2.0"
					<< std::endl
					<< " x_cells=10" << std::endl
					<< " y_cells=2" << std::endl
					<< " xmin=0.0" << std::endl
					<< " ymin=0.0" << std::endl
					<< " xmax=10.0" << std::endl
					<< " ymax=2.0" << std::endl
					<< " initial_timestep=0.04" << std::endl
					<< " timestep_rise=1.5" << std::endl
					<< " max_timestep=0.04" << std::endl
					<< " end_time=3.0" << std::endl
					<< " test_problem 1" << std::endl
					<< "*endclover" << std::endl;
			out_unit.close();
			g_in.open("clover.in");
		}
	}

	clover_barrier();

	if (parallel.boss) {
		g_out << std::endl
		      << "Initialising and generating" << std::endl
		      << std::endl;
	}

	read_input(g_in, parallel, config);

	clover_barrier();

//	globals.step = 0;
	config.number_of_chunks = parallel.max_task;


	auto globals = start(parallel, config, selectedDevice);

	clover_barrier(*globals);

	if (parallel.boss) {
		g_out << "Starting the calculation" << std::endl;
	}

	g_in.close();


	return globals;
}

