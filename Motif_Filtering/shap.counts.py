import argparse
import tensorflow
import kerasAC
import pandas as pd
import numpy as np
from scipy.special import softmax
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import *
from kerasAC.custom_losses import *
import math

import pdb
import argparse
import pickle
import tensorflow
from tensorflow.compat.v1.keras.backend import get_session
tensorflow.compat.v1.disable_v2_behavior()
import kerasAC
from kerasAC.generators.tiledb_predict_generator import *
from kerasAC.tiledb_config import *
from kerasAC.interpret.deepshap import *
from kerasAC.interpret.profile_shap import *
from kerasAC.helpers.transform_bpnet_io import *
from tensorflow.keras.models import load_model, model_from_json
#load the model!
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import *
from kerasAC.custom_losses import *

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--model_path")
    parser.add_argument("--onehot")
    parser.add_argument("--out_prefix")
    return parser.parse_args()

def get_interpretations(onehot, model, count_explainer):
    count_explanations=count_explainer.shap_values(onehot)[0]
    return count_explanations

def main():
    args = parse_args()
    model=model_from_json(open(args.model_path + '.arch', 'r').read())
    model.load_weights(args.model_path + '.weights')
    print("loaded_model")

    model_wrapper=(model.input, model.outputs[1])
    count_explainer=shap.DeepExplainer(model_wrapper,
                                       data=create_background_atac,
                                       combine_mult_and_diffref=combine_mult_and_diffref_1d)

    onehot = np.load(args.onehot)
    print("loaded onehot seq")

    batchsize=64

    final_shap = []

    for i in range(0, len(onehot), batchsize):
        print(str(int(i / batchsize)) + ' / ' + str(math.ceil(len(onehot)/batchsize)))
        batch = onehot[i:i+batchsize]
        print('Batch Shape: ' + str(batch.shape))

        shap_out = get_interpretations(batch, model, count_explainer)
        print('Output Shape: ' + str(shap_out.shape))

        if len(final_shap) == 0:
            final_shap = shap_out
        else:
            final_shap = np.concatenate((final_shap, shap_out))

    print('Final Shape: ' + str(final_shap.shape))
    np.save(args.out_prefix + '.count.shap.npy', final_shap)

if __name__=="__main__":
    main()

