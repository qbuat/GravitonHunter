import os
import shutil
for name in os.listdir("./"):
    if 'physics_Egamma' in name:
        # print name
        break_name = name.split(".")
        new_name = '.'.join([break_name[0], break_name[-1]])
        print name, '-->', new_name
        shutil.move(name, new_name)
