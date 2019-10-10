
from spack.cmd.install import *

def install_spec(cli_args, kwargs, abstract_spec, spec):
    """Do the actual installation."""

    # handle active environment, if any
    def install(spec, kwargs):
        env = ev.get_env(cli_args, 'install')
        if env:
            env.install(abstract_spec, spec, **kwargs)
            env.write()
        else:
            spec.package.do_install(**kwargs)

    try:
        if cli_args.things_to_install == 'dependencies':
            # Install dependencies as-if they were installed
            # for root (explicit=False in the DB)
            kwargs['explicit'] = False
            for s in spec.dependencies():
                install(s, kwargs)
            if create_host_config in dir(spec.package):
                spec.package.create_host_config(spec,spec.prefix)
        else:
            kwargs['explicit'] = True
            install(spec, kwargs)

    except spack.build_environment.InstallError as e:
        if cli_args.show_log_on_error:
            e.print_context()
            if not os.path.exists(e.pkg.build_log_path):
                tty.error("'spack install' created no log.")
            else:
                sys.stderr.write('Full build log:\n')
                with open(e.pkg.build_log_path) as log:
                    shutil.copyfileobj(log, sys.stderr)
        raise
