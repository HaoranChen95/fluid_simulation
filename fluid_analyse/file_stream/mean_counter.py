class mean_counter:
    def __init__(self):
        self.data = 0
        self.index = 0

    def __iadd__(self, data):
        self.index += 1
        self.data = (self.index - 1) / self.index * self.data + data / self.index
        return self
