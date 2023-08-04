import math
import random
import SM3
#生成n个随机数
def random_int_list(start, stop, length):
    start, stop = (int(start), int(stop)) if start <= stop else (int(stop), int(start))
    length = int(abs(length)) if length else 0
    return [random.randint(start, stop) for _ in range(length)]
 
 


#使用n个不同的数进行生日攻击
def brithday_attack(n):
    a=[]
    a=random_int_list(0,0xFFFFFFFF,n)
    print(a)
    b=[]
    for i in range(0,n-1):
        temp = SM3.Hash_sm3(str(a[i]),1)
        if temp not in b:
            b.append(temp)
        else:
            print("找到碰撞")
            
    
    

def print_res(i):
    
    print(i)
    print(SM3.Hash_sm3(str(i),1))
    print(SM3.Hash_sm3(str(i),1))


if __name__ == '__main__':

    #print(random_int_list(1, 100, 10))  
    brithday_attack(0xFFFFFFFF)

