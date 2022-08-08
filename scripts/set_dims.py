import argparse
import utils

#parse arguments supplied via command line
parser = argparse.ArgumentParser()
parser.add_argument('width', type=float, help='width in micrometers')
parser.add_argument('thickness', type=float, help='sample thickness in micrometers')
parser.add_argument('--gap-center', type=float, help='width in micrometers', default=9.0)
parser.add_argument('--params', type=str, help='name of parameter file', default="params.conf")
args = parser.parse_args()

#initialize information classes
geom = utils.Geometry(args.params)

left = str(geom.l_junc)
rght = str(geom.r_junc)
top =  str(geom.t_junc)
bot =  str(geom.b_junc)

with open(r'junc_template.eps', 'r') as file:
    data = file.read()
    data = data.replace('$LEFT', left)
    data = data.replace('$RGHT', rght)
    data = data.replace('$TOP', top)
    data = data.replace('$BOT', bot)
    print(data)
