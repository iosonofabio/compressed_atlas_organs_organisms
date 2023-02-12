# vim: fdm=indent
'''
author:     Fabio Zanini
date:       03/05/22
content:    Validate time points.
'''


def validate_correct_timepoint(text):
    '''Validate and correct a timepoint'''
    from text_recognition.assets import numbers

    text = text.strip(' ')
    if text.startswith('at '):
        text = text[3:]

    # The "e" in E18.5 seems hard to guess
    hard_to_guess_embryonic = False
    try:
        float(text)
        hard_to_guess_embryonic = True
    except ValueError:
        pass
    if hard_to_guess_embryonic:
        text = 'E'+text

    # Older ages are like 3m, 18m, etc.
    if text[0].isdigit():
        if text.endswith('months'):
            text = text[:-len('onths')]
        elif text.endswith('years'):
            text = text[:-len('ears')]

        if text[:-1].isdigit() and text[-1] in ('m', 'y'):
            return text
        else:
            return None

    # If not an older age, capitalize
    text = text[0].upper()+text[1:]

    # Not starting with a digit, so it shoul be e.g. E18.5/P1
    if not text[-1].isdigit():
        # Try to convert a string-style number (e.g. pseven)
        for ntext, ndigit in numbers[::-1]:
            if text.endswith(ntext):
                text = text[:-len(ntext)]+str(ndigit)
                break

    # If a dot, convert before the dot
    divider = text.find('.')
    if divider != -1:
        predot = text[1: divider]
        if not predot.isdigit():
            for ntext, ndigit in numbers[::-1]:
                if predot.endswith(ntext):
                    text = text[0] + str(ndigit) + text[divider:]
                    break
            else:
                return None

    try:
        float(text[1])
    except ValueError:
        return None

    return text
