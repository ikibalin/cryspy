"""The list of the user's packages (add_package) is kept in the user's
folder, as the package folder may be read-only."""
import os
import stat
import sys
import warnings

import pytest

import cryspy
from cryspy.H_functions_global import function_1_cryspy_objects as module


@pytest.fixture
def packages_files(tmp_path, monkeypatch):
    f_new = tmp_path / "user" / "packages.dat"
    f_old = tmp_path / "old" / "packages.dat"
    monkeypatch.setattr(module, "F_PACKAGES", str(f_new))
    monkeypatch.setattr(module, "F_OLD_PACKAGES", str(f_old))
    return f_new, f_old


def make_package(tmp_path):
    dir_package = tmp_path / "my_cryspy_package"
    dir_package.mkdir()
    (dir_package / "__init__.py").write_text(
        "def calc_nothing(a: float) -> float:\n"
        "    return a\n"
        "CRYSPY_CLASSES = []\n"
        "L_FUNCTION = [calc_nothing]\n")
    return str(dir_package)


def test_file_in_user_folder(monkeypatch):
    monkeypatch.delenv("CRYSPY_PACKAGES_FILE", raising=False)
    f_name = module.take_packages_file()
    assert not f_name.startswith(os.path.dirname(cryspy.__file__))
    assert os.path.basename(f_name) == "packages.dat"
    monkeypatch.setenv("CRYSPY_PACKAGES_FILE", "/some/file.dat")
    assert module.take_packages_file() == "/some/file.dat"


def test_add_load_delete(packages_files, tmp_path):
    f_new, _ = packages_files
    dir_package = make_package(tmp_path)
    assert cryspy.packages() == []
    cryspy.add_package(dir_package)
    cryspy.add_package(dir_package)  # once only
    assert [hh.strip() for hh in cryspy.packages()] == [dir_package]
    assert f_new.is_file()
    module.load_packages()
    assert "calc_nothing" in [f.__name__ for f in module.L_FUNCTION_ADD]
    module.delete_package(dir_package)
    assert cryspy.packages() == []


def test_old_file_is_read(packages_files, tmp_path):
    f_new, f_old = packages_files
    dir_package = make_package(tmp_path)
    f_old.parent.mkdir()
    f_old.write_text(dir_package + "\n" + str(tmp_path / "removed") + "\n")
    assert [hh.strip() for hh in cryspy.packages()] == [
        dir_package, str(tmp_path / "removed")]
    module.load_packages()
    # the folder that is not found is dropped; the list goes to the new file
    assert f_new.read_text().split() == [dir_package]


@pytest.mark.skipif(sys.platform == "win32" or os.geteuid() == 0,
                    reason="needs a read-only folder")
def test_read_only_folder_does_not_raise(packages_files, tmp_path):
    f_new, _ = packages_files
    f_new.parent.mkdir()
    f_new.write_text(str(tmp_path / "removed") + "\n")
    f_new.parent.chmod(stat.S_IRUSR | stat.S_IXUSR)
    f_new.chmod(stat.S_IRUSR)
    try:
        with warnings.catch_warnings(record=True) as l_warning:
            warnings.simplefilter("always")
            module.load_packages()
        assert any("not saved" in str(hh.message) for hh in l_warning)
    finally:
        f_new.parent.chmod(stat.S_IRWXU)
        f_new.chmod(stat.S_IRWXU)
