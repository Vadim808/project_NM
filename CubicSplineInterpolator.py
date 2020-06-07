import numpy as np

class CubicSplineInterpolator:
    
    def __init__(self,xGrid,fGrid):
        self.xGrid = xGrid #you will need it
        self.fGrid = fGrid #you will need it
        self.coeffs = self.ComputeCoefficients(xGrid,fGrid)
        
    def ComputeCoefficients(self, xGrid, fGrid):
        N = len(xGrid)
        h = np.append(np.append([1], xGrid[1:] - xGrid[:-1]), [1]) 
        I = np.append(np.append([0], fGrid[1:] - fGrid[:-1]), [0]) / h
        h[0], h[N] = 0, 0
        delta = [0] * (N + 1)
        delta[1] = -h[2] / (2 * (h[1] + h[2]))
        L = [0] * (N + 1)
        L[1] = 3/2 * (I[2] - I[1]) / (h[1] + h[2])
        for k in range(3, N + 1):
            delta[k - 1] = -h[k] / (2 * h[k - 1] + 2*h[k] + h[k - 1]*delta[k - 2])
            L[k - 1] = (3*I[k] - 3*I[k - 1] - h[k - 1] * L[k - 2]) / (2 * h[k - 1] + 2 * h[k] + h[k - 1] * delta[k - 2])
        c = [0] * (N + 1)
        for k in range(N, 1, -1):
            c[k-1] = delta[k - 1] * c[k] + L[k - 1]
        d = [0] * (N + 1)
        b = [0] * (N + 1)
        h = h + (h == 0)
        for k in range(1, N + 1):
            d[k] = (c[k] - c[k - 1])/ (3 * h[k]);
            b[k] = I[k] + (2 * c[k] * h[k] + h[k] * c[k-1]) / 3;
        return [b, c, d]
        
    def Compute(self, X):
        res = []
        for x in X:
            for i in range(1, len(self.xGrid)):
                if self.xGrid[i - 1] <= x <= self.xGrid[i]:
                    a = self.fGrid[i]
                    x_xk = x - self.xGrid[i]
                    b = self.coeffs[0][i]
                    c = self.coeffs[1][i]
                    d = self.coeffs[2][i]
                    res += [a + b * x_xk + c * x_xk ** 2 + d * x_xk ** 3]
                    break
        return res

    def Integrate(self, y, T=1):
        s = 0
        for i in range(1, len(self.xGrid)):
            if T < self.xGrid[i - 1]:
                break
            if self.xGrid[i] < y:
                continue
            if self.xGrid[i - 1] < y <= self.xGrid[i]:
                x0 = y - self.xGrid[i]
            if y <= self.xGrid[i - 1]:
                x0 = self.xGrid[i - 1] - self.xGrid[i]
            if self.xGrid[i - 1] < T <= self.xGrid[i]:
                x1 = T - self.xGrid[i - 1]
            else:
                x1 = 0
            a = self.fGrid[i]
            b = self.coeffs[0][i]
            c = self.coeffs[1][i]
            d = self.coeffs[2][i]
            A = a * x0 + b/2 * x0 ** 2 + c/3 * x0 ** 3 + d/4 * x0 ** 4
            B = a * x1 + b/2 * x1 ** 2 + c/3 * x1 ** 3 + d/4 * x1 ** 4
            s += B - A
        return s

    def Grad(self, X):
        res = []
        for x in X:
            for i in range(1, len(self.xGrid)):
                if self.xGrid[i - 1] <= x <= self.xGrid[i]:
                    a = self.fGrid[i]
                    x_xk = x - self.xGrid[i]
                    b = self.coeffs[0][i]
                    c = self.coeffs[1][i]
                    d = self.coeffs[2][i]
                    res += [b + 2 * c * x_xk + 3 * d * x_xk ** 2]
                    break
        return res

