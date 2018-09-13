

class Line:

    def __init__(self, vector1, vector2):
            self.p1 = vector1
            self.p2 = vector2

    def slope(self):
            '''Get the slope of a line segment'''
            (x1, y1), (x2, y2) = self.first, self.second
            try:
                return (float(y2)-y1)/(float(x2)-x1)
            except ZeroDivisionError:
                # The line is vertical
                return None

    def yintercept(self, slope):
            '''Get the y intercept of a line segment'''
            if slope != None:
                    x, y = self.first
                    return y - slope * x
            else:
                    return None

    def solve_for_y(self, x, slope, yintercept):
            '''Solve for Y cord using line equation'''
            if slope != None and yintercept != None:
                    return float(slope) * x + float(yintercept)
            else:
                    raise Exception('Can not solve on a vertical line')

    def solve_for_x(self, y, slope, yintercept):
            '''Solve for X cord using line equatio'''
            if slope != 0 and slope:
                    return float((y - float(yintercept))) / float(slope)
            else:
                    raise Exception('Can not solve on a horizontal line')