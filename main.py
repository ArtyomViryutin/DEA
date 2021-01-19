from dea.dea import dea

X = [
    [96, 16, 8],
    [110, 22, 1400],
    [100, 18, 1200],
    [125, 25, 1200],
    [120, 24, 1600]
]




Y = [
    [3800, 25, 8.0],
    [4600, 32, 8.5],
    [4400, 35, 8.0],
    [6500, 30, 10.0],
    [6000, 28, 9.0]
]



if __name__ == '__main__':
    results = dea(X, Y, RTS='vrs')
    print("Efficiency: ")
    for dmu, result in enumerate(results, start=1):
        print(f'DMU {dmu}: {result.fun}')
    # print(results[0])




