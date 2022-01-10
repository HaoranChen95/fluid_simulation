class color_map:
    def __init__(self):
        pass

    def line(self, start, end, x):
        resalt = []
        difference = []
        for x1, x2 in zip(start, end):
            difference.append(x2 - x1)
        for y0, y in zip(start[1:], difference[1:]):
            if difference[0]:
                resalt.append(y0 + (x - start[0]) * y / difference[0])
        return tuple(resalt)

    def blue_grey_red(self, x):
        if x < 0.5:
            return self.line(
                (0, 0x06 / 255, 0x9A / 255, 0xF3 / 255),
                (0.5, 0xD8 / 255, 0xD8 / 255, 0xD8 / 255),
                x,
            )
        else:
            return self.line(
                (0.5, 0xD8 / 255, 0xD8 / 255, 0xD8 / 255),
                (1, 0xFB / 255, 0x29 / 255, 0x43 / 255),
                x,
            )

    def blue_green_yellow_red(self, x):
        if x < 1 / 3:
            return self.line(
                (0, 0x06 / 255, 0x9A / 255, 0xF3 / 255),
                (1 / 3, 0xF1 / 255, 0xC4 / 255, 0x0F / 255),
                x,
            )
        elif x < 2 / 3:
            return self.line(
                (1 / 3, 0xF1 / 255, 0xC4 / 255, 0x0F / 255),
                (2 / 3, 0x27 / 255, 0xAE / 255, 0x60 / 255),
                x,
            )
        else:
            return self.line(
                (2 / 3, 0x27 / 255, 0xAE / 255, 0x60 / 255),
                (1, 0xFB / 255, 0x29 / 255, 0x43 / 255),
                x,
            )

    def color_str(self, rgb):
        r, g, b = [int(c * 255) for c in rgb]
        return f"#{r:02X}{g:02X}{b:02X}"
