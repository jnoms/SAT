from .struc_get_domains import visualize_PAE, parse_json_file
from .utils.misc import make_output_dir, talk_to_me


def plot_pae_main(args):
    talk_to_me("Parsing json file.")
    pae = parse_json_file(args.scores)[0]
    plt = visualize_PAE(pae)

    talk_to_me("Writing image")
    make_output_dir(args.out_image)
    plt.savefig(args.out_image)


if __name__ == "__main__":
    msg = "Call this script from sat.py, where there is argument parsing."
    raise ValueError(msg)
