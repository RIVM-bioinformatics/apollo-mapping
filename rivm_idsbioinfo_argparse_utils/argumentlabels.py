__doc__ = """ Functions that (temporarily) label the help_text of an argparse --flag with a special note to the user """

def label_expert_arg(add_arg_func):
    """ helper that wraps add_argument(), suffixing help to state this is an expert flag """
    def wrapper(*args, **kwargs):
        note = "EXPERT FLAG: use(full) only for pipeline developers"
        kwargs['help'] = f"{kwargs.get('help', '')} [{note}]"
        return add_arg_func(*args, **kwargs)
    return wrapper

def label_notimplemented_arg(add_arg_func):
    """ helper that wraps add_argument(), suffixing help to state this flag is not implemented yet """
    def wrapper(*args, **kwargs):
        note = "NotYetImplemented FLAG: program either ignores or will even fail when using this flag"
        kwargs['help'] = f"{kwargs.get('help', '')} [{note}]"
        return add_arg_func(*args, **kwargs)
    return wrapper

def label_deprecated_arg(add_arg_func):
    """ helper that wraps add_argument(), prefixing help to state this flag meanwhile became deprecated """
    def wrapper(*args, **kwargs):
        note = "DEPRECATED FLAG: program no longer supports it"
        kwargs['help'] = f"[{note}] {kwargs.get('help', '')}"
        return add_arg_func(*args, **kwargs)
    return wrapper