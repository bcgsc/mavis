import types


class Mock:
    def __init__(self, **kwargs):
        for attr, val in kwargs.items():
            setattr(self, attr, val)

    def bind_method(self, **kwargs):
        for attr, val in kwargs.items():
            val = types.MethodType(val, self)  # bind the method to self
            setattr(self, attr, val)

    def add_attr(self, attr, val):
        setattr(self, attr, val)


class MockFunction:
    def __init__(self, return_value):
        self.return_value = return_value

    def __call__(self, *pos, **kwargs):
        return self.return_value
