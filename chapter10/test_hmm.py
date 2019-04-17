from pathlib import Path
from stepik_hmm import fetch_hmm
import hmm


def test_viterbi():
    emissions, model = fetch_hmm(Path('test/testcase01.txt'))
    path = hmm.viterbi(emissions=emissions, model=model)
