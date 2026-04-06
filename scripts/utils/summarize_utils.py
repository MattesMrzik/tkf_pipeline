import re

def get_tree_params(path, tree_regex):
    match = re.search(tree_regex, path)
    return match.groupdict() if match else {}

def get_msa_sim_params(path, msa_sim_tools_config):
    for tool_name, tool_conf in msa_sim_tools_config.items():
        match = re.search(tool_conf["match_regex"], path)
        if match:
            params = {"msa_sim_tool": tool_name}
            params.update(match.groupdict())
            return params
    return {}

def get_inference_params(path, inference_tools_config):
    for inf_tool_name, inf_conf in inference_tools_config.items():
        match = re.search(inf_conf["match_regex"], path)
        if match:
            params = {"inference_tool": inf_tool_name}
            params.update(match.groupdict())
            return params
    return {}
