import sys
sys.path.append("train_utilities/")
import ann_utils
import json
import subprocess as sp
import os
import glob


class RegressionTrainExecutor:
    def __init__(self, mode, n_train=1000,n_test=500,n_valid=333):
        self.mode = mode
        self.n_train_examples = str(n_train)
        self.n_test_examples = str(n_test)
        self.n_valid_examples = str(n_valid)

    def cleanup(self):
        dirs_to_clean = ["./", "./train_data/", "./test_data/", "./validation_data/"]
        for i_dir in dirs_to_clean:
            path_to_files = glob.glob(os.path.join(i_dir, "*.txt"))
            for i_file in path_to_files:
                if os.path.isfile(i_file):
                    os.remove(i_file)

    def setup_files(self):
        ext = "txt"
        train_filename = "./train_data/{}_{}".format(self.mode, "TRAIN")
        valid_filename = "./validation_data/{}_{}".format(self.mode, "VALID")
        test_filename = "./test_data/{}_{}".format(self.mode, "TEST")
        comparison_filename= "./{}_comparison_file".format(self.mode)
        list_of_filenames = []
        for str_name in [train_filename, valid_filename, test_filename, comparison_filename]:
            list_of_filenames.append('.'.join([str_name, ext]))
        list_of_filenames = json.dumps(list_of_filenames)
        return list_of_filenames

    def setup_parameters(self):
        parameters_string = ''
        if self.mode.startswith("CE"):
            parameters = {
                    'kex' : (100.0,500.0),
                    'pb' : (0.05,0.15),
                    'deltaA' : (-10.0, -0.5),
                    'deltaB' : (0.5, 10.0),
                    'r1a' : (0.1, 4.0),
                    'r2a' : (5.0, 20.0),
                    'r2a' : (30.0, 200.0),
                    'tex' : (0.1, 0.4),
                    'b1field' : (20.0,80.0)
                    }
            parameters_string = json.dumps(parameters)
        elif self.mode.startswith("CP"):
            exchange_type = self.mode[4:-1]
            kex_ranges = (500.0, 1500.0) if exchange_type == "FAST" else (5.0, 600.0)
            parameters = {
                    'r2' : (5.0, 30.0),
                    'kex' : kex_ranges, # depends on exchange type
                    'dppmmin' : (1.0, 2.0),
                    'pa' : (0.8, 0.99),
                    'dppm' : (0.0, 3.0), # for slow exchange
                    'b0fields' : (50.0, 95.0),
                    }
            parameters_string = json.dumps(parameters)
        else:
            raise ValueError("No parameter values for chosen mode.\n")
        return parameters_string

    def main(self):
        file_names_json = self.setup_files()
        self.cleanup()
        parameters_json = self.setup_parameters()
        # Assuming the path for comdnmr is stored in PATH environment variable
        model_file = "train_utilities/cpmgTrain.py" if self.mode.startswith("CP") else "train_utilities/cestTrain.py"
        if os.path.isfile(model_file):
            sp.call(["comdnmr", model_file, file_names_json, parameters_json, self.n_train_examples, self.n_valid_examples, self.n_test_examples])
        elif os.path.isfile(model_file):
            sp.call(["comdnmr", model_file, file_names_json, parameters_json, self.n_train_examples, self.n_valid_examples, self.n_test_examples])
        else:
            raise LookupError("'{}' does not exist.\n".format(model_file))

#####################################################################################################
# args: mode, number of training examples, number of testing examples, number of validation examples

n_args = len(sys.argv)

mode = sys.argv[1].upper() if n_args >=  2 else '' # mode that corresponds to type of experiment 
train_examples = int(sys.argv[2]) if n_args > 2 and sys.argv[2].isdigit() else False
test_examples = int(sys.argv[3]) if n_args > 3 and sys.argv[3].isdigit() else False
valid_examples = int(sys.argv[4]) if n_args > 4 and sys.argv[4].isdigit() else False
available_modes = ["CEST", "CPMGFAST1", "CPMGSLOW1", "CPMGFAST2","CPMGSLOW2"]
assert mode in available_modes, "Mode is not valid. Options => '{}'\n".format(available_modes)

if mode and train_examples and test_examples and valid_examples:
    trainer = RegressionTrainExecutor(mode, train_examples, test_examples, valid_examples)
else:
    trainer = RegressionTrainExecutor(mode)

trainer.main()
#####################################################################################################




