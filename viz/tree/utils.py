import re

def get_tree_params(path, tree_regex):
    match = re.search(tree_regex, path)
    return match.groupdict() if match else {}
