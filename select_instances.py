def name_from_path(path):
    words = path.split("/")
    fullname = words[-1]
    return fullname.split(".")[0]

def instance_selection():
    easy_instances = []
    remaining_instances = []
    too_hard = []
    with open("singular_results.out", "r") as f:
        for line in f.readlines():
            if "timed out" in line:
                name = name_from_path(line.split()[0])
                too_hard.append(name)
                continue
            else:
                words = line.split()
                path = words[0]
                time = float(words[1])
                name = name_from_path(path)
                if time < 0.1:
                    easy_instances.append(name)
                elif time <= 2:
                    remaining_instances.append(name)
                else:
                    too_hard.append(name)
    return easy_instances, remaining_instances, too_hard
