import json
import os

MODEL_INF_DIR = "model_inference"

def load_params_json(path):
    if not os.path.exists(path):
        return {}
    with open(path, "r") as f:
        data = json.load(f)
    return data

