from colour import Color


def dynamic_label_color(color):
    """
    calculates the luminance of a color and determines if a black or white label will be more contrasting
    """
    f = Color(color)
    if f.get_luminance() < 0.5:
        return '#FFFFFF'
    else:
        return '#000000'


class LabelMapping:
    def __init__(self, **kwargs):
        self._mapping = dict()
        self._reverse_mapping = dict()
        for k, v in kwargs.items():
            self[k] = v

    def __setitem__(self, key, value):
        if key in self._mapping:
            raise KeyError('duplicate key: keys must be unique', key)
        if value in self._reverse_mapping:
            raise KeyError('duplicate value: values must be unique', value)
        self._mapping[key] = value
        self._reverse_mapping[value] = key

    def items(self):
        return self._mapping.items()

    def __getitem__(self, key):
        return self._mapping[key]

    def __len__(self):
        return len(self._mapping.keys())

    def get_key(self, value):
        return self._reverse_mapping[value]

    def set_key(self, key, value):
        if key in self._mapping:
            current_value = self._mapping[key]
            if value == current_value:
                return
            elif value in self._reverse_mapping:
                raise KeyError('duplicate value: values must be unique', value)
            del self._mapping[key]
            del self._reverse_mapping[current_value]
        elif value in self._reverse_mapping:
            raise KeyError('duplicate value: values must be unique', value)
        self[key] = value

    def add(self, value, prefix=''):
        if value in self._reverse_mapping:
            return self._reverse_mapping[value]
        i = 1
        while True:
            key = '{}{}'.format(prefix, i)
            if key not in self._mapping:
                self[key] = value
                break
            i += 1
        return self._reverse_mapping[value]
