def convert_time(t):
    if t < 60:
        t_unit = 's'
    elif t < 3600:
        t /= 60
        t_unit = 'm'
    else:
        t /= 3600
        t_unit = 'h'

    return t, t_unit
