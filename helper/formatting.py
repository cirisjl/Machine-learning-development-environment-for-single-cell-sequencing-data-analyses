def convert_arr_to_sql_str(arr):
    if type(arr[0]) == str:
        arr = ["'" + str(x) + "'" for x in arr]
    else:
        arr = [str(x) for x in arr]
    return "(" + ", ".join(arr) + ")"
