def overlap(start_1, end_1, start_2, end_2):
    result = min([end_1, end_2]) - max([start_1, start_2]) + 1
    if result >= 0:
        return result
    return 0
