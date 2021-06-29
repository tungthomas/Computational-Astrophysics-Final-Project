import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# parameters need to specified first
figure_size = 1
speed = 500

count = 0
with open('data.txt', 'r') as file:
    line = file.readlines()
    for i in line:
        if i.find('---------') == -1:
            count += 1
        else:
            break
        # print(i)
# print(line[0])

np = count-1
steps = int(len(line)/(np+2))
pos = [[[0 for k in range(3)] for j in range(np)] for i in range(steps)]
time = []

for i in range(steps):
    for j in range(np+2):
        n = (np+2)*i + j
        if(j == 0):
            time.append(float(line[n][5:]))
        elif(j == np+1):
            continue
        else:
            sp = line[n].find('.')
            a = float(line[n][sp+7:sp+15])
            b = float(line[n][sp+17:sp+25])
            c = float(line[n][sp+27:sp+35])
            # print(a, b, c, (i, j))
            pos[i][j-1][0] = a
            pos[i][j-1][1] = b
            pos[i][j-1][2] = c
        # print(line[n])

print('number of particle is ', np)
print(steps)

dt = time[1]-time[0]
end_time = dt*steps

# plot the figure

nstep_per_image = 1


# create figure
fig = plt.figure(figsize=(6, 6), dpi=100)
ax = plt.axes(xlim=(-figure_size, figure_size),
              ylim=(-figure_size, figure_size))

balls = []
for i in range(np):
    b, = ax.plot([], [], 'ro', ms=5)
    balls.append(b)

text = ax.text(0.0, 24, '', fontsize=16, color='black',
               ha='center', va='center')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_aspect('equal')
ax.tick_params(top=True, right=True, labeltop=True, labelright=True)
# ax.add_artist(plt.Circle((0.0, 0.0), r, color='b', fill=False))

count = 0


def init():
    for i in range(np):
        balls[i].set_data(pos[0][i][0], pos[0][i][1])
    text.set(text='')
    return balls, text


def update(i):
    global count
    count += 1
    print(count)
    for i in range(np):
        balls[i].set_data(pos[count][i][0], pos[count][i][1])
    text.set(text='hello world')
    return balls, text


# create movie
nframe = steps-1
anim = animation.FuncAnimation(fig, func=update, init_func=init,
                               frames=nframe, interval=speed, repeat=False)
plt.show()
