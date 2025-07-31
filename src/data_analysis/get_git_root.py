import git


def get_git_root(path):
    git_repo = git.Repo(path, search_parent_directories=True)
    return git_repo.working_dir
