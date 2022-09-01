import argparse
import utils

LEFT_STR = "$LEFT"
RGHT_STR = "$RGHT"
TOP_STR = "$TOP"
BOT_STR = "$BOT"
DELTA_TOP_STR = "$DELTA_TOP"
DELTA_BOT_STR = "$DELTA_BOT"
FIELD_AMP_STR = "$FIELD_AMP"
LEFT_LEN = len(LEFT_STR)
RGHT_LEN = len(RGHT_STR)
TOP_LEN = len(TOP_STR)
BOT_LEN = len(BOT_STR)
DELTA_TOP_LEN = len(DELTA_TOP_STR)
DELTA_BOT_LEN = len(DELTA_BOT_STR)
FIELD_AMP_LEN = len(FIELD_AMP_STR)

#parse arguments supplied via command line
parser = argparse.ArgumentParser()
parser.add_argument('width', type=float, help='width in micrometers')
parser.add_argument('thickness', type=float, help='sample thickness in micrometers')
parser.add_argument('--in-file', type=str, help='file to read as input', default='junc_template.geom')
parser.add_argument('--gap-center', type=float, help='width in micrometers', default=9.0)
parser.add_argument('--field-amp', type=float, help='amplitude of the field in V.A^-1', default=1.7)
parser.add_argument('--params', type=str, help='name of parameter file', default="params.conf")
args = parser.parse_args()

def read_word(s, i):
    '''Read the word from the string s starting at index i
    '''
    if i >= len(s):
        return '', 0
    #find the first non-whitespace character
    j = i
    while s[j] == ' ' or s[j] == '\t':
        j += 1
    #hacky way to treat operators as words
    if s[j] == '-' or s[j] == '+' or s[j] == '*' or s[j] == '/':
        return s[j:j+1], j+1
    #iterate until we find a word terminator
    while j < len(s):
        if s[j] == ' ' or s[j] == '\t' or s[j] == '(' or s[j] == ')' or s[j] == '[' or s[j] == ']' or s[j] == ',' or s[j] == '-' or s[j] == '+' or s[j] == '*' or s[j] == '/':
            return s[i:j], j
        j += 1
    return s[i:], j

#initialize information classes
geom = utils.Geometry(args.params, args.width, args.thickness)

'''left = str(geom.l_junc)
rght = str(geom.r_junc)
top =  str(geom.t_junc)
bot =  str(geom.b_junc)
delta_top = str(geom.t_junc - 1.0)
delta_bot = str(geom.t_junc + 1.0)
field_amp_meep = str(geom.mks_e_field_to_meep_field(args.field_amp*10e10))'''
left = geom.l_junc
rght = geom.r_junc
top =  geom.t_junc
bot =  geom.b_junc
delta_top = geom.t_junc - 1.0
delta_bot = geom.t_junc + 1.0
field_amp_meep = geom.mks_e_field_to_meep_field(args.field_amp*10e10)

def parse_value(s, tok, i):
    '''parse the string s[i:] and the matching token tok to find a real value. i should be the index of the first character after the token
    '''
    val = 0.0
    try:
        val = float(tok)
    except ValueError:
        pass

    j = i
    if tok == LEFT_STR:
        val = left
        j = i + LEFT_LEN
    elif tok == RGHT_STR:
        val = rght
        j = i + RGHT_LEN
    elif tok == DELTA_TOP_STR:
        val = delta_top
        j = i + DELTA_TOP_LEN
    elif tok == DELTA_BOT_STR:
        val = delta_bot
        j = i + DELTA_BOT_LEN
    elif tok == TOP_STR:
        val = top
        j = i + TOP_LEN
    elif tok == BOT_STR:
        val = bot
        j = i + BOT_LEN
    elif tok == FIELD_AMP_STR:
        val = field_amp_meep
        j = i + FIELD_AMP_LEN
    #read the following word
    op = ''
    next_word, j = read_word(s, j)
    tv = 0.0
    if len(next_word) == 1:
        op = next_word
        next_word, j = read_word(s, j)
        tv, tj = parse_value(s, next_word, j)
    #this doesn't obey order of operations
    if op == '-':
        val = val - tv
        j = tj
    elif op == '+':
        val = val + tv
        j = tj
    elif op == '*':
        val = val * tv
        j = tj
    elif op == '/':
        val = val / tv
        j = tj
    return val, j

with open(args.in_file, 'r') as fl:
    data = fl.read()
    out_data = ''
    last_var_ind = 0
    var_ind = data.find('$')
    while var_ind >= 0:
        var_name, skip = read_word(data, var_ind)
        val, skip = parse_value(data, var_name, var_ind)
        out_data = out_data + data[last_var_ind:var_ind] + str(val)
        #look for the next variable and the end
        last_var_ind = skip
        var_ind = data.find('$', last_var_ind)
    #we have to write everything that comes after the last variable
    out_data = out_data + data[last_var_ind:]
    '''data = data.replace('$LEFT', left)
    data = data.replace('$RGHT', rght)
    data = data.replace('$DELTA_TOP', delta_top)
    data = data.replace('$DELTA_BOT', delta_bot)
    data = data.replace('$TOP', top)
    data = data.replace('$BOT', bot)
    data = data.replace('$FIELD_AMP', field_amp_meep)'''
    print(out_data)
