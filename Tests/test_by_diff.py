import argparse, os, re

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description='Unit test for runs with fixed seed, which compares base results to dynamically generated ones by absolute diff')
    parser.add_argument('--run_directory', '-r', help='directory cd should be cded to upon running', required=False, default = os.getcwd() + "/Examples/")
    parser.add_argument('--traitrelax_program_path', '-t', help='full path to ',
                        required=False, default = os.getcwd() + "/TraitRELAX/TraitRELAX")
    parser.add_argument('--parameter_file_path', '-p', help='path to input parameter file for the program',
                        required=True)
    parser.add_argument('--output_directory', '-o', help='directory that holds the base files and the output files that should be compared to them', required=False, default="")
    parser.add_argument('--comparison_filenames', '-c', help='file names for comparison', required=False, default=["optimization_output.res", "parameter_estimates.res"])

    args = parser.parse_args()
    run_directory = args.run_directory
    traitrelax_program_path = args.traitrelax_program_path
    parameter_file_path = args.parameter_file_path
    example_name_regex = re.compile("([^\/]*?)\.bpp", re.DOTALL)
    example_name = example_name_regex.search(parameter_file_path).group(1)
    output_directory = args.output_directory
    if output_directory == "":
        output_directory = run_directory + "data/" + example_name + "/"
    print("example_name: ", example_name)
    comparison_filenames = args.comparison_filenames
    if type(comparison_filenames) != list:
        comparison_filenames = comparison_filenames.split(",")

    # execute the program
    res = os.chdir(run_directory)
    res = os.system(traitrelax_program_path + " param=" + parameter_file_path)

    # make sure that the output files were created
    for filename in comparison_filenames:
        if not os.path.exists(output_directory + filename):
            print("file " + output_directory + filename + " does not exist")
            exit(1)

    # make sure the base files exist
    for filename in comparison_filenames:
        if not os.path.exists(output_directory + "base_" + filename):
            print("file " + output_directory +  "base_" + filename + " does not exist")
            exit(1)

    # compare base files to output files
    for filename in comparison_filenames:
        output_file_path =  output_directory + filename
        base_file_path = output_directory +  "base_" + filename
        with open(output_file_path, "r") as output_file:
            output_content = output_file.read()
        with open(base_file_path, "r") as base_file:
            base_content = base_file.read()
        if output_content != base_content:
            print("base and output fies for " + filename + " are different")
            res=os.system("diff " + base_file_path + " " + output_file_path)
            exit(1)

    # remove output files
    for filename in comparison_filenames:
        res = os.system("rm -r " + output_directory + filename)
