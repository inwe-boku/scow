from pathlib import Path

from scow import paths


def test_path_helpers_create_missing_directories(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "repo", tmp_path)

    results = paths.results_dir()
    figures = paths.figures_dir()
    opt = paths.opt_dir()

    assert results == tmp_path / "data" / "results"
    assert figures == tmp_path / "figures"
    assert opt == tmp_path / "opt"
    assert results.is_dir()
    assert figures.is_dir()
    assert opt.is_dir()
    assert (tmp_path / "data").is_dir()


def test_result_file_helpers_write_under_created_results_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "repo", tmp_path)

    expected = tmp_path / "data" / "results"
    assert paths.coefs_path("2014", "full").parent == expected
    assert paths.loglik_path("2014", "full").parent == expected
    assert paths.soco_path("2014").parent == expected
