from catplot.asimov import _is_approved


def test_latest_review_wins():
    assert _is_approved({"review": [{"status": "APPROVED"}]})
    assert not _is_approved({"review": [{"status": "APPROVED"}, {"status": "REJECTED"}]})
    assert _is_approved({"review": [{"status": "REJECTED"}, {"status": "APPROVED"}]})
    assert not _is_approved({"review": []})
    assert not _is_approved({})
