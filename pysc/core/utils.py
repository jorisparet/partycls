import re
from pysc.particle import aliases

def _standardize_condition(condition):
    """
    Check that the condition is correctly formated (i.e <attr> _operator_ <val>).

    Parameters
    ----------
    condition : str
        condition.

    Raises
    ------
    ValueError
        If condition is not valid or if the <attr> is not recognized).

    Returns
    -------
    condition : str
        A standardized condition.

    """
    regexp = re.search('(\w+\.?\w+)\s?(<|<=|==|>=|>)\s?([\'|\"]?\w+[\'|\"]?)', condition)
    if regexp:
        attr = regexp.group(1)
        operator = regexp.group(2)
        value = regexp.group(3)
        # if attribute is an alias, replace it with full attribute
        if not(attr.startswith('particle')):
            if attr in aliases:
                attr = aliases[attr]
            else:
                raise ValueError('attribute "{}" is not recognized'.format(attr))
        # reformat condition
        condition = '{} {} {}'.format(attr, operator, value)
        return condition
    else:
        raise ValueError('"{}" is not a valid condition'.format(condition))