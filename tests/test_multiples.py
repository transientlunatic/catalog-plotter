import numpy as np

from catplot.multiples import bounded_kde, jensen_shannon_divergence, split_pages


def test_reflected_kde_integrates_to_one_on_bounded_support():
    rng = np.random.default_rng(1)
    samples = rng.uniform(0, 1, 5000)
    grid = np.linspace(0, 1, 501)
    density = bounded_kde(samples, grid, bounds=(0, 1))
    assert abs(np.trapz(density, grid) - 1) < 0.02
    # No depletion at the edges of a flat distribution
    assert density[0] > 0.8 and density[-1] > 0.8


def test_kde_has_no_support_outside_bounds():
    rng = np.random.default_rng(2)
    samples = np.abs(rng.normal(size=3000))
    grid = np.linspace(-1, 4, 200)
    density = bounded_kde(samples, grid, bounds=(0, None))
    assert np.all(density[grid < 0] == 0)


def test_jsd_limits():
    x = np.linspace(-10, 10, 2000)
    p = np.exp(-0.5 * (x + 6) ** 2)
    q = np.exp(-0.5 * (x - 6) ** 2)
    assert jensen_shannon_divergence(p, p) == 0
    assert abs(jensen_shannon_divergence(p, q) - 1) < 1e-6
    assert jensen_shannon_divergence(p, q) == jensen_shannon_divergence(q, p)


def test_split_pages_is_balanced():
    pages = split_pages(list(range(28)), 14)
    assert [len(p) for p in pages] == [14, 14]
    pages = split_pages(list(range(29)), 14)
    assert sorted(len(p) for p in pages) == [9, 10, 10]
    assert split_pages(list(range(5)), 0) == [list(range(5))]
