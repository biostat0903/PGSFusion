import random

males = open('male.id', 'r').readlines()
females = open('female.id', 'r').readlines()

males = [x.strip() for x in males]
females = [x.strip() for x in females]

validation1_males = random.sample(males, 25000)
validation1_females = random.sample(females, 25000)

print(len(males))

with open('validatioin1.id', 'w') as f:
    for x in validation1_males:
        f.write(x + ' ' + x + '\n')
    for x in validation1_females:
        f.write(x + ' ' + x + '\n')

males = list(set(males).difference(set(validation1_males)))
females = list(set(females).difference(set(validation1_females)))

validation2_males = random.sample(males, 25000)
validation2_females = random.sample(females, 25000)

print(len(males))

with open('validatioin2.id', 'w') as f:
    for x in validation2_males:
        f.write(x + ' ' + x + '\n')
    for x in validation2_females:
        f.write(x + ' ' + x + '\n')

males = list(set(males).difference(set(validation2_males)))
females = list(set(females).difference(set(validation2_females)))

test_males = random.sample(males, 25000)
test_females = random.sample(females, 25000)

print(len(males))

with open('test.id', 'w') as f:
    for x in test_males:
        f.write(x + ' ' + x + '\n')
    for x in test_females:
        f.write(x + ' ' + x + '\n')

males = list(set(males).difference(set(test_males)))
females = list(set(females).difference(set(test_females)))

ref_males = random.sample(males, 1000)
ref_females = random.sample(females, 1000)

print(len(males))

with open('ref.id', 'w') as f:
    for x in ref_males:
        f.write(x + ' ' + x + '\n')
    for x in ref_females:
        f.write(x + ' ' + x + '\n')

males = list(set(males).difference(set(ref_males)))
print(len(males))
