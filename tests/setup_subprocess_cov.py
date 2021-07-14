import os
import sys

for p in sys.path:
    if p.endswith('site-packages'):
        pth_file = os.path.join(p, 'subprocess-coverage.pth')
        print('writing path file:', pth_file)
        with open(pth_file, 'w') as fh:
            fh.write('import coverage\n\ncoverage.process_startup()\n')
        break
