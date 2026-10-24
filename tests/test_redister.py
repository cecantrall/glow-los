from glow_los.rt import redister


def test_redister_runtime_is_temporary_and_independent_of_retrieval_package(monkeypatch, tmp_path):
    monkeypatch.setattr(redister.tempfile, "gettempdir", lambda: str(tmp_path))
    with redister._redister_runtime_dir() as runtime:
        assert runtime.parent == tmp_path / "glow-los"
        assert runtime.is_dir()
        (runtime / "redister.out").write_text("test")
    assert not runtime.exists()
