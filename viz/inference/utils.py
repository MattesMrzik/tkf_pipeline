import re

def get_inference_params(path, inference_tools_config):
    for inf_tool_name, inf_conf in inference_tools_config.items():
        match = re.search(inf_conf["match_regex"], path)
        if match:
            params = {"inference_tool": inf_tool_name}
            params.update(match.groupdict())
            return params
    return {}

