def cutter(string:str):
    """
    This function cuts a string in multiple parts.
    """
    cutted_string = string.split("/")
    cutted_string = [cutted.strip() for cutted in cutted_string]
    return cutted_string