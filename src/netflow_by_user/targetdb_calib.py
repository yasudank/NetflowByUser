import os
import sys
from sshtunnel import SSHTunnelForwarder
from sqlalchemy import create_engine, text, sql
import toml
from astropy.table import Table
import pandas as pd
from tqdm import tqdm

def get_search_radius(
    fp_radius_degree = 260.0 * 10.2 / 3600,  # "Radius" of PFS FoV in degree
    fp_fudge_factor = 1.5  # fudge factor for search widths
):
    # search radius 
    search_radius = fp_radius_degree * fp_fudge_factor
    print("search_radius is %f degree." % search_radius)

    return search_radius

def get_config(config_fn):
    with open(config_fn, "r") as f:
        config = toml.load(f)
    print(config)

    return config

def get_ppcList(config):
    fn = os.path.join(config["input"]["dir"], config["input"]["fn_ppcList"]) 
    ppcList = Table.read(fn)
    print("There are %d pointings." % len(ppcList))
    print("ppcList read from %s" % fn)

    return ppcList

def get_centerList(config):
    fn = os.path.join(config["input"]["dir"], config["input"]["fn_ppcList"]) 
    ppcList = Table.read(fn)
    print("There are %d pointings." % len(ppcList))
    print("ppcList read from %s" % fn)

    return ppcList

def get_sky(config, engine, ra, dec, search_radius, frac=0.1):
    tablename = 'sky'
    version = config['targetdb']['sky']['version']

    query_string = f"""SELECT *
    FROM {tablename}
    WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {search_radius}) 
          AND (version = '{version}') 
          AND (random() < {frac})
    """
    conn = engine.connect()
    query = conn.execute(sql.text(query_string))
    df = pd.DataFrame(query.fetchall())
    conn.close()

    return df

def get_fluxstd(config, engine, ra, dec, search_radius, frac=0.2):
    tablename = 'fluxstd'
    version = config['targetdb']['fluxstd']['version']
    min_prob_f_star = config['targetdb']['fluxstd']['min_prob_f_star']
    input_catalog_id = config['targetdb']['fluxstd']['input_catalog_id']

    query_string = f"""SELECT *
    FROM {tablename}
    WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {search_radius}) 
          AND (prob_f_star BETWEEN {min_prob_f_star} AND 1.0
          OR is_fstar_gaia = True)
          AND (version = '{version}')
          AND (random() < {frac});
    """
    
    conn = engine.connect()
    query = conn.execute(sql.text(query_string))
    df = pd.DataFrame(query.fetchall())
    conn.close()

    return df

def main(engine_b):

    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    else:
        config_file = 'config_targetdb_cosmos.toml'

    # set seach radius in degree
    search_radius = get_search_radius()

    # read config from config_file
    config = get_config(config_file)

    # read ppcList
    #ppcList = get_ppcList(config)
    #ppc_code_list, ra_list, dec_list = ppcList['ppc_code'], ppcList['ppc_ra'], ppcList['ppc_dec']

    # read center list
    ppcList = get_centerList(config)
    print(ppcList)
    ppc_code_list, ra_list, dec_list = ppcList['ppc_code'], ppcList['ppc_ra'], ppcList['ppc_dec']
    #sys.exit(1)

    sky_dir = os.path.join(config['output']['dir'], "sky")
    if not os.path.exists(sky_dir):
        os.system(f'mkdir -p {sky_dir}')

    fluxstd_dir = os.path.join(config['output']['dir'], "fluxstd")
    if not os.path.exists(fluxstd_dir):
        os.system(f'mkdir -p {fluxstd_dir}')
    
    for ppc_code, ra, dec in tqdm(
        zip(ppc_code_list, ra_list, dec_list),
        total=len(ppc_code_list),
        desc="Processing pointings"
    ):
        
        df = get_sky(config, engine_b, ra, dec, search_radius, 1.0)
        outfn = os.path.join(sky_dir, f"{ppc_code}.ecsv")

        table = Table.from_pandas(df)
        table.write(outfn, format="ascii.ecsv", overwrite=True)

        tqdm.write("%s: %d sky selected." % (ppc_code, len(df)))
        tqdm.write('write to %s'%outfn)

        df = get_fluxstd(config, engine_b, ra, dec, search_radius, 1.0)
        outfn = os.path.join(fluxstd_dir, f"{ppc_code}.ecsv")

        table = Table.from_pandas(df)
        table.write(outfn, format="ascii.ecsv", overwrite=True)

        tqdm.write("%s: %d fluxstd selected." % (ppc_code, len(df)))
        tqdm.write('write to %s'%outfn)
        
# 接続設定ファイルのロード
db_config_file = os.environ.get(
    "DB_CONFIG_FILE",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 
        "db_config.toml"
    )
)

if not os.path.exists(db_config_file):
    print(f"Error: Database config file not found at '{db_config_file}'")
    sys.exit(1)

with open(db_config_file, "r") as f:
    db_conf = toml.load(f)

# 1. SSH接続情報 (hostA)
# ssh-agent を起動し、鍵を登録しておく
# eval "$(ssh-agent -s)"
# ssh-add /home/yasuda/.ssh/xxxxxx
SSH_HOST = db_conf["ssh"]["host"]
SSH_USER = db_conf["ssh"]["user"]
SSH_PKEY_PATH = db_conf["ssh"]["pkey_path"]

# 2. PostgreSQL接続情報 (hostB)
DB_B_HOST = db_conf["db_b"]["host"]
DB_B_PORT = db_conf["db_b"]["port"]
DB_B_USER = db_conf["db_b"]["user"]
DB_B_PASSWORD = db_conf["db_b"]["password"]
DB_B_NAME = db_conf["db_b"]["name"]

# SSHトンネルの確立
with SSHTunnelForwarder(
    (SSH_HOST, 22),
    ssh_username=SSH_USER,
    # ssh_pkeyを指定しないことで、自動的にssh-agentの鍵が使用されます
    remote_bind_addresses=[
        (DB_B_HOST, DB_B_PORT)
    ]
) as server:
    
    # 各接続先に対応するローカルポートを取得 (tunnel_bindingsの値は (local_ip, local_port) のタプルです)
    local_port_b = server.tunnel_bindings[(DB_B_HOST, DB_B_PORT)][1]
    
    print(f"hostB トンネル完了 -> localhost:{local_port_b}")

    # 4. hostB への接続エンジン作成
    db_url_b = f"postgresql+psycopg2://{DB_B_USER}:{DB_B_PASSWORD}@127.0.0.1:{local_port_b}/{DB_B_NAME}"
    engine_b = create_engine(db_url_b)
    
    main(engine_b)

    engine_b.dispose()

# withブロックを抜けると、SSHトンネルは自動的にクローズされます
print("SSHトンネルをクローズしました。")
