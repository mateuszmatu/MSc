dicts = {0: 'A', 0: 'B'}
keys = [10, 12]
values = ["A", "B"]
print(dicts)
for i in range(len(keys)):
    dicts[keys[i]] = values[i]
print(dicts)