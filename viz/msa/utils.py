import re

def get_msa_sim_params(path, msa_sim_tools_config):
    for tool_name, tool_conf in msa_sim_tools_config.items():
        match = re.search(tool_conf["match_regex"], path)
        if match:
            params = {"msa_sim_tool": tool_name}
            params.update(match.groupdict())
            return params
    return {}
