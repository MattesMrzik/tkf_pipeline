import json
import os

from viz.utils import get_tool_type_dir

MODEL_INF_DIR = get_tool_type_dir("param_inf")

def load_params_json(path):
    if not os.path.exists(path):
        return {}
    with open(path, "r") as f:
        data = json.load(f)
    return data

