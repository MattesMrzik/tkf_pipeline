import numpy as np

def prob_stop_links_is_success(lam, mu):
    return 1.0 - (lam / mu)

def prob_of_zero_length(lam, mu):
    return prob_stop_links_is_success(lam, mu)

def expected_length(lam, mu, r):
    return (lam / mu) / ((1 - (lam / mu)) * (1 - r))

def simulate_tkf_root_length(n_experiments, lam, mu, r):
    stop_links_is_success = prob_stop_links_is_success(lam, mu)
    prob_stop_fragment_is_success = 1.0 - r

    results = []

    for _ in range(n_experiments):
        num_links = np.random.geometric(stop_links_is_success) - 1
        
        total_length = 0
        for _ in range(num_links):
            # Equivalent to sample_fragment_length
            # Fragments in TKF are at least 1 unit long
            frag_len = np.random.geometric(prob_stop_fragment_is_success)
            total_length += frag_len
            
        results.append(total_length)
    return np.array(results)

def plot_hist_vs_theory(results, lam, mu, r):
    import matplotlib.pyplot as plt
    counts = np.bincount(results)
    probs_emp = counts / counts.sum()

    k = np.arange(len(probs_emp))

    # theoretical geometric
    success =  (1 - lam / mu) * (1-r)
    probs_theory = lam/mu * (1 - success) ** (k - 1) * success

    # plot
    plt.figure()
    plt.bar(k, probs_emp, alpha=0.5, label="Simulation")
    plt.plot(k[1:], probs_theory[1:], label="Theory")

    prob_of_zero_length = prob_stop_links_is_success(lam, mu)
    plt.scatter(0,  prob_of_zero_length, color='red', label="Zero length prob (theory)")

    zero_length_fraction = np.mean(results == 0)
    plt.scatter(0, zero_length_fraction, color='blue', label="Zero length prob (simulation)")

    expectation = expected_length(lam, mu, r)
    simulated_mean = np.mean(results)

    plt.xlabel("Total length")
    plt.ylabel("Probability")
    plt.legend()
    # the expected prob for having length zero is having zero num_links which is 
    plt.title("TKF root length: simulation vs analytical; expected mean = {:.2f}, simulated mean = {:.2f}".format(expectation, simulated_mean))
    plt.show()

def sim_and_plot(lam, mu, r, n_experiments=1000000):
    results = simulate_tkf_root_length(n_experiments, lam, mu, r)
    plot_hist_vs_theory(results, lam, mu, r)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Simulate TKF root length and compare to theory.")
    parser.add_argument("-l", type=float, required=True, help="Lambda parameter for TKF model")
    parser.add_argument("-m", type=float, required=True, help="Mu parameter for TKF model")
    parser.add_argument("-r", type=float, required=True, help="R parameter for TKF model")
    parser.add_argument("-n", type=int, help="Number of simulations to run")
    args = parser.parse_args()

    if args.n:
        print("Prob of zero length: {:.4f}".format(prob_of_zero_length(args.l, args.m)))
        print("Expected length: {:.4f}".format(expected_length(args.l, args.m, args.r)))
        sim_and_plot(args.l, args.m, args.r, args.n)
    else:
        print("Prob of zero length: {:.4f}".format(prob_of_zero_length(args.l, args.m)))
        print("Expected length: {:.4f}".format(expected_length(args.l, args.m, args.r)))
