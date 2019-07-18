# purpose of this script is to generate random dummy data to help train and test an artificial neural network
# that is to receive a field(s) along with R2(eff) point data as input and output the values for the dispersion curve that 
# describes CPMG Fast or Slow exchange. 

# Note: This script is being called from runner.sh.

from org.comdnmr.fit.calc import CPMGEquation
import random as rd
import sys
sys.path.append("./train_utilities")
sys.path.append("./train_utilities/java_ann_files/NeuralNet.jar")
#sys.path.append("/Users/teddycolon/NetBeansProjects/NeuralNet/dist/NeuralNet.jar")
#sys.path.append("/Users/teddycolon/NetBeansProjects/NeuralNet/dist/lib/ojalgo-47.2.0.jar")
from ANN import NeuralNetworkUtils as netutils
import ann_utils
import os
import subprocess as sp
import json

# There are 2 CPMG modes to generate a data for:
# 1_ CPMG FAST EXCHANGE: pars => R2, KEX, dPPMmin
# 2_ CPMG SLOW EXCHANGE: pars => R2, KEX, REX, pA, dPPM
################################################################################################
# GLOBAL SCOPE
# all the command line args below:
FILENAMES_JSON = json.loads(sys.argv[1])
PARAMETERS_JSON = json.loads(sys.argv[2])
EXPERIMENT_TYPE = FILENAMES_JSON[0].split('/')[-1].split('.')[0].split('_')[0]
NEXAMPLESLIST = [int(val) for val in sys.argv[3:] if val.isdigit()] # train, valid, test
R2EFF_MAX_VALUE = 150.0
INFOLIST = [] # list containing all the info before writing to file

################################################################################################

# Switcher to grab parameter labels based on type of experiment 
def get_params(cpmg_type):
    switcher = {
        'CPMGFAST' : ("kex", "r2", "dppmmin"),
        'CPMGSLOW' : ("kex", "pa", "r2", "dppm")
    }
    return switcher.get(cpmg_type, None)

# Random parameter generator (uniform distribution)
def generate_random_values(params):
    randomized_params = []
    for param in params:
         low, high = PARAMETERS_JSON.get(param)
         random_value = rd.uniform(low, high)
         randomized_params.append(random_value)
    return randomized_params

# Random field generator (uniform distribution)
def generate_random_fields(field_max, field_min, n_fields):
    border_offset = 1.0 # offset value to prevent generated field value from being the minimum or maximum value.
    offset_field_min = field_min + border_offset
    initial_field = float(rd.randint(offset_field_min, 80.0))
    current_fields = []
    if n_fields >= 1:
        offset_bt_fields = 10.0
        offset_field_min = initial_field + offset_bt_fields
        offset_field_max = field_max - border_offset
        current_fields = [initial_field if n == 0 else float(rd.randint(offset_field_min, offset_field_max)) for n in range(n_fields)]
    return current_fields

# Scale the randomly generated parameter values between 0-1
def scale_pars(parameters_to_scale, param_labels):
    scaled_pars = []
    for i, param in enumerate(param_labels):
        low, high = PARAMETERS_JSON[param]
        scaled_pars.append(ann_utils.scaleValue(parameters_to_scale[i], high, low))
    return scaled_pars

# Calculate R2,eff values for one or multiple fields depending on experiment
def calculate_ann_input(data_block):
    global R2EFF_MAX_VALUE # included due to possible changes in SCALE_VALUE
    id_num = 0
    x_arr = data_block.get("x_vals")
    pars = data_block.get("par")
    fields = data_block.get("fields")
    n_fields = data_block.get("n_fields")
    mode = data_block.get("mode")
    y_values = []

    if mode.startswith("FA"):
        map_par = CPMGEquation.CPMGFAST.makeMap(1)
        for j in range(n_fields):
            for val in x_arr:
                output_val = CPMGEquation.CPMGFAST.calculate(pars, map_par[0], [val], id_num, fields[j])
                y_values.append(output_val)
    elif mode.startswith("SL"):
        map_par = CPMGEquation.CPMGSLOW.makeMap(1)
        for j in range(n_fields):
            for val in x_arr:
                output_val = CPMGEquation.CPMGSLOW.calculate(pars, map_par[0], [val], id_num, fields[j])
                y_values.append(output_val)

    # check if there's a value higher than the scale value to report it at the end.
    #print y_values
    if max(y_values) > R2EFF_MAX_VALUE:
        R2EFF_MAX_VALUE = max(y_values)

    return y_values

# Bruce guesser 
def current_guesser(exchange_type, x_values, y_values, field):
    """
    - exchange_type (String) : determines whether to execute the fast or slow exchange guesser code
    - x_values (array<floats>)
    - y_values (array<floats>)
    - map_arr (2D array<integers>)
    - id_nums (array<integers>)
    - n_id (integer)
    - field (float)
    """
    from org.comdnmr.fit.calc import CalcRDisp
    from org.comdnmr.fit.calc import DataUtil
    from org.comdnmr.fit.calc import CoMDPreferences
    import math

    n_id = 0
    id_nums = [n_id]*len(y_values)

    if exchange_type == "FAST":
        map_arr = CPMGEquation.CPMGFAST.makeMap(1)
        n_pars = CalcRDisp.getNPars(map_arr)
        guesses = [0.0]*n_pars
        kex_sum = 0.0

        id_example_numex = 0
        while (id_example_numex < len(map_arr)):

            min_y = DataUtil.getMinValue(y_values, id_nums, id_example_numex)
            max_y = DataUtil.getMaxValue(y_values, id_nums, id_example_numex)
            mean = DataUtil.getMeanValue(y_values, id_nums, id_example_numex)
            v_mid = DataUtil.getMidValue(y_values, x_values, id_nums, id_example_numex)
            r2 = min_y * 0.95
            rex = max_y - min_y
            if rex < 0.0:
                rex = 0.0
            guesses[map_arr[id_example_numex][1]] = r2
            tau_mid = 1.0 / (2.0 * v_mid)
            kex = 1.915 / (0.5 * tau_mid)
            dPPMMinRad = math.sqrt(4.0 * rex / (field * field) * kex)
            dPPMMin = dPPMMinRad / (2.0 * math.pi)
            guesses[map_arr[id_example_numex][2]] = dPPMMin
            if rex >= 0.0:
                kex_sum += kex
            
            id_example_numex += 1

        guesses[0] = kex_sum / len(map_arr)
        if guesses[0] > CoMDPreferences.getCPMGMaxFreq():
            guesses[0] = CoMDPreferences.getCPMGMaxFreq() * 0.9

        return guesses
    elif exchange_type == "SLOW":
        map_arr = CPMGEquation.CPMGSLOW.makeMap(1)
        n_pars = CalcRDisp.getNPars(map_arr)
        guesses = [0.0] * n_pars
        kex_sum = 0.0
        pa = 0.95

        id_example_numex = 0
        while (id_example_numex < len(map_arr)):
            min_y = DataUtil.getMinValue(y_values, id_nums, id_example_numex)
            max_y = DataUtil.getMaxValue(y_values, id_nums, id_example_numex)
            mean = DataUtil.getMeanValue(y_values, id_nums, id_example_numex)
            v_mid = DataUtil.getMidValue(y_values, x_values, id_nums, id_example_numex)
            r2 = min_y * 0.95
            rex = max_y - r2
            tau_mid = 1.0 / (2.0 * v_mid)
            kex = 1.915 / (0.5 * tau_mid)
            if kex > CoMDPreferences.getCPMGMaxFreq():
                kex = CoMDPreferences.getCPMGMaxFreq() * 0.9

            dw2 = rex / (pa * (1.0 - pa)) * kex
            dPPM = math.sqrt(dw2) / (2.0 * math.pi) / field
            guesses[map_arr[id_example_numex][2]] = r2
            guesses[map_arr[id_example_numex][3]] = dPPM
            kex_sum += kex

            id_example_numex += 1

        guesses[0] = kex_sum / len(map_arr)
        guesses[1] = pa
        return guesses

    else:
        raise ValueError("The exchange type '{}' is not valid.".format(exchange_type))

# Package the example lines prepared and stored in INFOLIST
def write_example_lines(write_file):
    for line in INFOLIST:
        write_file.write(line)

# Train neural network 
def train_network(train_file, validation_file, scale_values, save_ann_directory, in_out_neurons):
    from ANN import CPMGTrain

    if os.path.isfile(train_file) and os.path.isfile(validation_file):
        print "-- training the network!"
        train_obj = CPMGTrain(train_file, validation_file)
        train_obj.setScales(scale_values)
        n_input_neurons, n_output_neurons = in_out_neurons

        train_label = EXPERIMENT_TYPE # label for the network being trained
        result = train_obj.cpmgANNRun(train_label, n_input_neurons, n_output_neurons)
        if result:
            trained_ann = result.getFirst()
            ann_info = result.getSecond()

            saved_ann_file = netutils.saveNeuralNetwork(save_ann_directory, trained_ann, ann_info)
            print "-- network saved!"
            return saved_ann_file
        else:
            raise ValueError("Result object after training is not what expected.")
    else:
        raise ValueError("train and/or validation file do not exist.")

# Test metwork
def test_network(saved_ann_file, ann_input):
    print "-- testing network!"
    return  netutils.getSavedNetworkOutput(saved_ann_file, ann_input)

# Helper function to call trainer and tester
def performance_eval(data_block):
    print "-- performing evaluation!"

    input_for_ann = data_block.get("ann_input")
    target_output = data_block.get('ann_target')
    saved_ann_file = data_block.get("saved_ann")
    if saved_ann_file:
        ann_guess = test_network(saved_ann_file, input_for_ann)
    else:
        raise ValueError("Saved ANN file string is empty")
    print("-- test complete!")

    revert_scales = data_block.get("reverter")
    label_list = data_block.get("label_list")

    # for current
    x_vals = data_block.get("xvals")
    y_vals = data_block.get("yvals")
    exc_type = data_block.get("mode")
    curr_guess = []
    curr_guess_per_field = []
    fields = data_block.get("fields")
    len_of_fields = len(fields)
    for field in fields:
        curr_guess = current_guesser(exc_type, x_vals, y_vals,field)
        if len_of_fields > 1:
            curr_guess_per_field.append(curr_guess)

    if len_of_fields > 1:
        curr_guess = ann_utils.getMean(curr_guess_per_field)

    # scale current guesser output
    scaled_curr_guess = []
    for i, label in enumerate(label_list):
        low, high = PARAMETERS_JSON.get(label)
        scaled_curr_guess.append(ann_utils.scaleValue(curr_guess[i], high, low))


    compare_file_obj = data_block.get("comp_file_obj")

    # scaled:
    compare_file_obj.write("-- scaled information: \n")
    ann_utils.comparisonFile(label_list, target_output, scaled_curr_guess, ann_guess, compare_file_obj)

    # unscaled:
    r_target_output = [0.0]*len(target_output)
    r_ann_guess = [0.0]*len(ann_guess)
    for i, label in enumerate(label_list):
        high, low = revert_scales[label]
        r_target_output[i] = ann_utils.revertScale(target_output[i], high, low)
        r_ann_guess[i] = ann_utils.revertScale(ann_guess[i], high, low)

    compare_file_obj.write("-- unscaled information: \n")
    ann_utils.comparisonFile(label_list, r_target_output, curr_guess, r_ann_guess, compare_file_obj)
        
    ann_rmse = ann_utils.getRMSE(target_output, ann_guess)
    curr_rmse = ann_utils.getRMSE(target_output, scaled_curr_guess)
    rmse_list = [curr_rmse, ann_rmse]
    compare_rmse_string = "curr_rmse : %f, ann_rmse : %f\n" % (curr_rmse, ann_rmse)
    compare_rmse_string += "******************************************************\n"
    compare_file_obj.write(compare_rmse_string)
    return rmse_list

# Main function
def main_func():
    ann_file_objects  = [ann_utils.makeFile(str(str_file)) for str_file in FILENAMES_JSON[:-1]]
    comparison_file_object = ann_utils.makeFile(FILENAMES_JSON[-1])
    n_fields = 1 if ann_file_objects[0].name.__contains__("1") else 2
    field_min, field_max = PARAMETERS_JSON.get('b0fields')
    r2eff_high, r2eff_low = 150.0, 0.0
    x_vals = [10, 20, 50, 100, 200, 400, 600, 800, 1000, 1100]
    is_test = False
    experiment = EXPERIMENT_TYPE[:-1]
    is_fast = True if experiment.__contains__("FAST") else False
    mode = "FAST" if is_fast else "SLOW"
    parameter_labels = get_params(experiment)
    saved_ann_dir = "./saved_networks"
    for i_file, n_examples in enumerate(NEXAMPLESLIST):
        if ann_file_objects[i_file].name.__contains__("TEST"):
            # check if test
            is_test = True
            rmse_list = []
            list_of_scales = []
            revert_scale_dict = {}
            for label in parameter_labels:
                low, high = PARAMETERS_JSON.get(label)
                revert_scale_dict[label] = (high, low)
                boundaries = ','.join([str(high),str(low)])
                labeled_bounds = ';'.join([label, boundaries])
                list_of_scales.append(labeled_bounds)
            list_of_scales.append("r2eff;{0},{1}".format(r2eff_high, r2eff_low))
            list_of_scales.append("fields;{0},{1}".format(field_max, field_min))
            scale_vals = ':'.join(list_of_scales) 
            train_file_path = FILENAMES_JSON[0]
            validation_file_path = FILENAMES_JSON[1]
            in_out_neurons = [11] if n_fields == 1 else [22] # 11 || 22 input neurons depending on n_fields
            in_out_neurons.append(3) if is_fast else in_out_neurons.append(4) # 3 || 4 depending on exchange type
            saved_ann_file = train_network(train_file_path, validation_file_path, scale_vals, saved_ann_dir, in_out_neurons)
            
        current_fields = generate_random_fields(field_max, field_min, n_fields)
        i_example = 1
        list_of_fields = [current_fields]
        while i_example <= n_examples:
            target_rand_params = generate_random_values(parameter_labels)
            data_block = {'x_vals' : x_vals, 'par' : target_rand_params, 'fields' : current_fields, 'n_fields' : n_fields, 'mode' : mode}
            y_vals = calculate_ann_input(data_block)
            if max(y_vals) > r2eff_high:
                continue
            else:
                i_example += 1
            scaled_input_fields = [(ann_utils.scaleValue(field, field_max, field_min)) for field in current_fields]
            original_y_vals = y_vals
            scaled_target_pars = scale_pars(target_rand_params, parameter_labels)

            line = ''
            if not is_test:
                y_vals = ann_utils.generateNoise(0, r2eff_high*0.007, y_vals)
            scaled_y_vals = [(ann_utils.scaleValue(r2eff_value, r2eff_high, r2eff_low)) for r2eff_value in y_vals]
            ann_input_vals = scaled_input_fields + scaled_y_vals
            line += ann_utils.generateOutputString(ann_input_vals, scaled_target_pars)
            line += '\n'
            INFOLIST.append(line)

            if is_test:
                ann_data_block = {
                    "saved_ann" : saved_ann_file,
                    "xvals" : x_vals,
                    "yvals" : original_y_vals,
                    "label_list" : parameter_labels,
                    "ann_input" : ann_input_vals,
                    "ann_target" : scaled_target_pars,
                    "comp_file_obj" : comparison_file_object,
                    "fields" : current_fields,
                    "reverter" : revert_scale_dict,
                    "mode" : mode
                    }
                comparison_file_object.write("\t\tTest number ({})\n".format(i_example- 1))
                rmse = performance_eval(ann_data_block)
                rmse_list.append(rmse)
                list_of_fields.append(current_fields)
                current_fields = generate_random_fields(field_max, field_min, n_fields)
                print "\t\t(_Test number %d_)" % (i_example - 1)
                print "field : %d" % (current_fields[0])
                print "\torig_R2eff\tscaled_R2eff"
                for ind, _ in enumerate(original_y_vals):
                    print "\t%f\t%f" % (original_y_vals[ind], scaled_y_vals[ind])
                print "\t-----------------------------------------------"
                print "\ttarget_param\tscaled_target_param:"
                for ind, _ in enumerate(target_rand_params):
                    print "\t%f\t%f" % (target_rand_params[ind], scaled_target_pars[ind])
            list_of_fields.append(current_fields)
            current_fields = generate_random_fields(field_max, field_min, n_fields)
        if is_test and rmse_list and list_of_fields:
            ann_utils.closeFile(comparison_file_object)
            # rmse_pair : curr_rmse, ann_rmse
            str_rmse = ','.join([':'.join([str(rmse_pair[0]), str(rmse_pair[1])]) for rmse_pair in rmse_list])
            # save plot data
            with open("{}_PLOTDATA.txt".format(EXPERIMENT_TYPE), 'w') as f_plot:
                bar_plot_input_string = "{} {} {}".format(mode, str(n_fields), str_rmse)
                f_plot.write(bar_plot_input_string)
            sp.Popen(["python", "train_utilities/displayBarPlot.py", mode, str(n_fields), str_rmse])



        write_example_lines(ann_file_objects[i_file])
        print("File created : %s" % (ann_file_objects[i_file].name))
        ann_utils.closeFile(ann_file_objects[i_file])

#########

main_func()
#print("Largest R2,eff value = {}".format(R2EFF_MAX_VALUE))
