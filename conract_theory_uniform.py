from scipy import optimize
# linspace의 데이터 추출 시각화
import matplotlib.pyplot as plt
import math
import numpy as np

def border() :
    print()
    print("=================================================")
    print()

def sigma(K, N, value):
    return sum(value[n] for n in range(K, N))

def positive_number(x) :
    return x if x > 0 else None

def incentive_X_delay_graph(Contract_item = []) :
    graph_color = ['c', 'b', 'g', 'r', 'm', 'y', 'k']
    for i in range(len(Contract_item)) :
        x = list(filter(positive_number, Contract_item[i].Delay))
        y = list(filter(positive_number, Contract_item[i].Incentive))
        label = "r=" + str(Contract_item[i].Gamma) + ",r'=" + str(Contract_item[i].Gamma_Prime)
        plt.plot(x,y, graph_color[i], label = label)
    plt.axis([0, 6, 0, 3.5])
    plt.xlabel("Latency")
    plt.ylabel("Incentive$")
    plt.title("Relationships between delay and incentive")
    plt.grid(True)
    plt.legend()
    # plt.savefig("contract1.png", dpi=350)
    plt.show()


def HubTpye_X_HubUtility_graph(Hub_type={}, N = 0, st=0, dt=0) :
    graph_color = ['c', 'b', 'g', 'r', 'm', 'y', 'k']
    for i in Hub_type :
        if i < st :
           continue
        elif i > dt :
            break
        x = list(Hub_type.keys())
        y = Hub_type[i]
        label = "type=" + str(i)
        color = graph_color.pop(0)
        plt.scatter(x,y, s = 40, c = color, alpha=0.5)
        plt.plot(x,y, color, label = label)
    plt.axis([10, N, -3, 6])
    plt.xlabel("Type Hub node")
    plt.ylabel("Utility Hub node")
    plt.title("Relationships between Type and Utility (Hub Node)")
    plt.grid(True)
    plt.legend()
    # plt.savefig("contract2.png", dpi=350)
    plt.show()



class Contract_item :
    def __init__(self, name):
        self.name = name
        self.Incentive = [0]
        self.Client_U = 0     # total client Utility
        self.Client_U_I = [0]
        self.Hub_U = [0]
        self.Hub_type_U = {}
        self.Delay = [0]
        self.Delay_Inverse = [0]
        self.Omega = [0]        # delay addition function
        self.Theta = []
        self.P = []
        self.Gamma = 0          # Incentive wegiht parameter for client
        self.Minus = -1
        self.N = 0              # total Theta number
        self.Gamma_Prime = 0    # unit resource cost for hub nodes

    def set_Theta_number(self, n):
        self.N = n
        self.Theta = [i for i in range(n+1)]
        self.P = [1/n for _ in range(n+1)]

        self.print_N()
        self.print_Theta()
        self.print_P()

    def set_Gamma(self, gamma):
        self.Gamma = gamma

        self.print_Gamma()

    def set_Gamma_Prime(self, gamma_prime):
        self.Gamma_Prime = gamma_prime

        self.print_Gamma_Prime()

    def print_N(self):
        print("{}'s N value : {}".format(self.name, self.N))

    def print_Theta(self):
        print("{}'s Theta value : {}".format(self.name, self.Theta))

    def print_P(self):
        print("{}'s P value : {}".format(self.name, self.P))

    def print_Gamma(self):
        print("{}'s Gamma value :       {}".format(self.name, self.Gamma))

    def print_Gamma_Prime(self):
        print("{}'s Gamma_Prime value : {}".format(self.name, self.Gamma_Prime))

    def state_print(self, name, value):
        for i in range(len(value)):
            print("{}'s {} index {} : {}".format(self.name, i, name, value[i]))

    def set_client_utility_i(self, i):
        def optimizeQ(x):
            current_P = self.P[i]
            V_theta_q = self.Theta[i] * math.log(1 + x)
            P_sigma = sigma(i + 1, self.N + 1, self.P)
            if i != self.N:
                Lambda = self.Theta[i] * math.log(1 + x) - self.Theta[i + 1] * math.log(1 + x)
            else:
                Lambda = 0

            return (self.Minus * ((current_P / self.Gamma_Prime * V_theta_q) +
                             (Lambda / self.Gamma_Prime * P_sigma) -
                             (current_P * self.Gamma * x)))

        return optimizeQ

    def set_omega(self):
        for k in range(1, self.N+1):
            if k == 1:
                self.Omega.append(0)
            else:
                self.Omega.append(self.Theta[k] * (math.log(1 + self.Incentive[k]) - math.log(1 + self.Incentive[k - 1])))

    def set_delay(self, i):
        delay = (self.Theta[1] * math.log(1 + self.Incentive[1]) + sigma(1, i + 1, self.Omega)) / self.Gamma_Prime
        self.Delay_Inverse.append(delay)
        if delay != 0:
            delay = 1 / delay
        return (delay)

    def set_client_Utility(self):
        return sum(self.P[i] * (self.Delay_Inverse[i] - self.Gamma * self.Incentive[i]) for i in range(1, self.N + 1))

    def set_hub_Utility(self):
        for i in range (1, self.N+1) :
            self.Hub_U.append(self.Theta[i] * math.log(1 + self.Incentive[i]) - self.Gamma_Prime * self.Delay_Inverse[i])

    def set_hub_type_Utility(self):
        for i in range (1, self.N+1) :
            self.Hub_type_U[self.Theta[i]] = []
            for j in range (1, self.N+1) :
                self.Hub_type_U[self.Theta[i]].append(self.Theta[i] * math.log(1 + self.Incentive[j]) - self.Gamma_Prime * self.Delay_Inverse[j])

    def execute(self,x0,bnds):
        # 1_ for
        for i in range(1, self.N+1) :
            temp = self.set_client_utility_i(i)
            result = optimize.minimize(temp, x0, method="TNC", bounds=bnds, options={'maxiter': 1000})
            self.Client_U_I.append(result.fun)
            self.Incentive.append(result.x)

        self.state_print("client utility", self.Client_U_I)
        border()
        self.state_print("incentive", self.Incentive)
        border()
        self.set_omega()
        self.state_print("Omega", self.Omega)
        border()

        #2_ for
        for i in range(1, self.N+1) :
            self.Delay.append(self.set_delay(i))

        self.state_print("Delay", self.Delay)
        border()

        self.set_hub_Utility()
        self.state_print("Hub Utility", self.Hub_U)
        border()

        self.set_hub_type_Utility()
        # border()

        print("**=====================================================")
        print()
        self.Client_U = self.set_client_Utility()
        print("Client_U : {}".format(self.Client_U))
        print()
        print("**=====================================================")



contract_item1 = Contract_item("item1")
contract_item2 = Contract_item("item2")
contract_item3 = Contract_item("item3")

contract_item1.set_Theta_number(20)
contract_item2.set_Theta_number(20)
contract_item3.set_Theta_number(20)
border()
contract_item1.set_Gamma(1)
contract_item1.set_Gamma_Prime(5)
contract_item2.set_Gamma(1.5)
contract_item2.set_Gamma_Prime(5)
contract_item3.set_Gamma(1)
contract_item3.set_Gamma_Prime(4.6)
border()

x0 = [0.1]          # init
bnds = [(0.0,40.0)] # bound

contract_item1.execute(x0, bnds)
contract_item2.execute(x0, bnds)
contract_item3.execute(x0, bnds)

incentive_X_delay_graph([contract_item1,contract_item2, contract_item3])
HubTpye_X_HubUtility_graph(contract_item3.Hub_type_U,N = contract_item3.N, st = 14, dt = 19)

