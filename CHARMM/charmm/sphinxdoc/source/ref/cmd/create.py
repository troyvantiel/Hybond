import os
import re

template = []
with open('command.rst.template', 'r') as t:
    for line in t.readlines():
        template.append(line)

with open('list', 'r') as f:
    for command in f.readlines():
        upper = command.strip().upper()
        lower = upper.lower()
        output = []
        for line in template:
            tmp = re.sub(r'REPLACEME', upper, line)
            tmp = re.sub(r'replaceme', lower, tmp)
            output.append(tmp)
        with open('%s.rst'%lower, 'w') as out:
            out.write("".join(output))

