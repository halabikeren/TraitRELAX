import argparse, re, os

def parse_inference(path):
    null_logl_regex = re.compile(
        "Null model fitting.*Overall Log likelihood\.*\:\s*(-\d*\.?\d*).*?\**Alternative model fitting",
        re.MULTILINE | re.DOTALL)
    alternative_logl_regex = re.compile(
        "Alternative model fitting.*Overall Log likelihood\.*\:\s*(-\d*\.?\d*).*?\**Performing statistical test",
        re.MULTILINE | re.DOTALL)
    with open(path, "r") as infile:
        content = infile.read()
    log_of_likelihood_ratio = float("nan")
    try:
        null_logl = float(null_logl_regex.search(content).group(1))
        alternative_logl = float(alternative_logl_regex.search(content).group(1))
        log_of_likelihood_ratio = alternative_logl - null_logl
    except Exception as e:
        print("couldn't parse " + path + " due to error: ", e)
    return  log_of_likelihood_ratio


if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(
        description="performs a statistical test given likelihood of a model under tha null and alternative "
                    "hypotheses based on an empirical distribution of likelihood ratios")
    parser.add_argument('--background_inference_dir', '-i',
                        help='directory that holds the inference output for datasets simulated under the null '
                             'hypothesis based on the parameter inference of TraitRELAX',
                        required=True)
    parser.add_argument('--null_log_likelihood', '-nll',
                        help='the value of null log likelihood obtained by TraitRELAX program', required=True)
    parser.add_argument('--alternative_log_likelihood', '-all',
                        help='the value of alternative log likelihood obtained by TraitRELAX program', required=True)
    args = parser.parse_args()
    background_inference_dir = args.background_inference_dir
    null_log_likelihood = float(args.null_log_likelihood)
    alternative_log_likelihood = float(args.alternative_log_likelihood)

    # parse the empirical inferences to obtain a distribution of log likelihood values
    log_of_likelihood_ratios = []
    for path in os.listdir(background_inference_dir):
        log_of_likelihood_ratios.append(parse_inference(background_inference_dir+path))
    log_of_likelihood_ratios = list(filter((float("nan")).__ne__, log_of_likelihood_ratios))

    # get the placement of the likelihood ratio of interest across the distribution
    original_log_likelihood_ratio = alternative_log_likelihood - null_log_likelihood
    log_of_likelihood_ratios.append(original_log_likelihood_ratio)
    log_of_likelihood_ratios.sort(reverse=True)
    original_placement = log_of_likelihood_ratios.index(original_log_likelihood_ratio)
    pvalue = original_placement / (len(log_of_likelihood_ratios)-1) # if the placement is below the 5% placement across the distribution, reject the null
    if pvalue <= 0.05:
        print("The null hypothesis has been rejected with a pvale of ", pvalue)
    else:
        print("The null hypothesis wasn't rejected. P-value is ", pvalue)



