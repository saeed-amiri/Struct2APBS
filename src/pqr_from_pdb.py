"""
make the data from the pdb and itp
"""

import sys
import pandas as pd

import logger
import itp_to_df
import pdb_to_df
import parse_charmm_data


class PdbToPqr:
    """
    preapre the file with positions and charges and radii
    """

    def __init__(self,
                 log: logger.logging.Logger
                 ) -> None:
        self.initiate(log)

    def initiate(self,
                 log: logger.logging.Logger
                 ) -> pd.DataFrame:
        """get all the infos"""
        itp: pd.DataFrame = itp_to_df.Itp('APT_COR.itp').atoms
        pdb: pd.DataFrame = pdb_to_df.Pdb(sys.argv[1], log).pdb_df
        charmm: pd.DataFrame = \
            parse_charmm_data.ParseData('CHARMM.DAT', log).radius_df
        pdb_with_charges: pd.DataFrame = self.get_charges(pdb, itp)
        pdb_with_charge_radii: pd.DataFrame = \
            self.set_radii(pdb_with_charges, charmm, itp)
        pdb_df: pd.DataFrame = self.add_chain_identifier(pdb_with_charge_radii)
        pqr_df: pd.DataFrame = self.mk_pqr_df(pdb_df)
        self.write_pqr(pqr_df)

    def get_charges(self,
                    pdb: pd.DataFrame,
                    itp: pd.DataFrame
                    ) -> pd.DataFrame:
        """get charges of the atoms in the pdb file"""
        aptes_with_charges: pd.DataFrame = self.set_aptes_charges(pdb, itp)
        cores_with_charges: pd.DataFrame = self.set_cores_charges(pdb, itp)
        return pd.concat(
            [cores_with_charges, aptes_with_charges], axis=0, ignore_index=True
            )

    def set_radii(self,
                  pdb_with_charges: pd.DataFrame,
                  charmm: pd.DataFrame,
                  itp: pd.DataFrame
                  ) -> pd.DataFrame:
        """set the radii for all"""
        aptes_df: pd.DataFrame = \
            self.set_radii_for_aptes(pdb_with_charges, charmm)
        cores_df: pd.DataFrame = \
            self.set_radii_for_cores(pdb_with_charges, charmm, itp)
        return pd.concat([cores_df, aptes_df], axis=0, ignore_index=True)

    @staticmethod
    def add_chain_identifier(pdb_df: pd.DataFrame
                             ) -> pd.DataFrame:
        """add the column"""
        chain_identifier_map = {
            'COR': 'A',
            'APT': 'B'
        }
        pdb_df['chain_id'] = \
            pdb_df['residue_name'].map(chain_identifier_map)
        return pdb_df

    @staticmethod
    def mk_pqr_df(pdb_with_charge_radii: pd.DataFrame
                  ) -> pd.DataFrame:
        """prepare df in the format of the pqr file"""
        columns: list[str] = ['records',
                              'atom_id',
                              'atom_name',
                              'residue_name',
                              'chain_id',
                              'residue_number',
                              'x', 'y', 'z', 'charge', 'radius']
        float_columns: list[str] = ['x', 'y', 'z', 'charge', 'radius']
        df_i: pd.DataFrame = pdb_with_charge_radii[columns].copy()
        df_i[float_columns] = df_i[float_columns].astype(float)
        return df_i

    @staticmethod
    def write_pqr(pqr_df: pd.DataFrame
                  ) -> None:
        """writing the pqr to a file"""
        with open('APT_COR.pqr', 'w', encoding='utf8') as f_w:
            for _, row in pqr_df.iterrows():
                line = f"ATOM  {row['atom_id']:>5} " \
                       f"{row['atom_name']:<4} " \
                       f"{row['residue_name']:<3} " \
                       f"{row['chain_id']:>1} " \
                       f"{row['residue_number']:>5} " \
                       f"{row['x']:>8.3f}" \
                       f"{row['y']:>8.3f}" \
                       f"{row['z']:>8.3f} " \
                       f"{row['charge']:>7.4f} " \
                       f"{row['radius']:>6.4f}\n"
                f_w.write(line)
            f_w.write('TER\n')
            f_w.write('END\n')

    @staticmethod
    def set_radii_for_aptes(pdb_with_charges: pd.DataFrame,
                            charmm: pd.DataFrame
                            ) -> pd.DataFrame:
        """get the radius from the file"""
        aptes_radii: pd.DataFrame = charmm[charmm['resname'] == 'APT'].copy()
        aptes_df = pd.DataFrame = \
            pdb_with_charges[pdb_with_charges['residue_name'] == 'APT'].copy()
        aptes_df = pd.merge(aptes_df,
                            aptes_radii[['atom_name', 'radius']],
                            on='atom_name',
                            how='left')
        return aptes_df

    @staticmethod
    def set_radii_for_cores(pdb_with_charges: pd.DataFrame,
                            charmm: pd.DataFrame,
                            itp: pd.DataFrame
                            ) -> pd.DataFrame:
        """get the radius from the file"""
        cores_radii: pd.DataFrame = charmm[charmm['resname'] == 'COR'].copy()
        cores_df = pd.DataFrame = \
            pdb_with_charges[pdb_with_charges['residue_name'] == 'COR'].copy()
        cores_df['atomtype'] = itp[itp['resname'] == 'COR']['atomtype'].copy()
        cores_df = pd.merge(cores_df,
                            cores_radii[['atom_name', 'radius']],
                            left_on='atomtype',
                            right_on='atom_name',
                            how='left')
        cores_df['atom_name'] = cores_df['atom_name_x'].copy()
        return cores_df.drop(
            ['atomtype', 'atom_name_x', 'atom_name_y'], axis=1)

    def set_aptes_charges(self,
                          pdb: pd.DataFrame,
                          itp: pd.DataFrame
                          ) -> pd.DataFrame:
        """get charges for different aptes longs"""
        itp_df: pd.DataFrame = itp[itp['resname'] == 'APT']
        apt_df: pd.DataFrame = pdb[pdb['residue_name'] == 'APT']
        pro_charges: pd.DataFrame
        unpro_charges: pd.DataFrame
        pro_charges, unpro_charges = self.get_aptes_charges(itp_df)
        dfs_with_charges: list[pd.DataFrame] = []
        for _, group in apt_df.groupby('residue_number'):
            if len(group) == 13:
                dfs_with_charges.append(self.add_charge_column_to_df(
                    group, pro_charges))
            elif len(group) == 12:
                dfs_with_charges.append(self.add_charge_column_to_df(
                    group, unpro_charges))
        return pd.concat(dfs_with_charges, axis=0, ignore_index=True)

    def set_cores_charges(self,
                          pdb: pd.DataFrame,
                          itp: pd.DataFrame
                          ) -> pd.DataFrame:
        """set the charges for the core atoms of the nano particles
        since only the charges on the shell is important, only those
        are checked
        """
        itp_df: pd.DataFrame = itp[itp['resname'] == 'COR'].copy()
        cor_df: pd.DataFrame = pdb[pdb['residue_name'] == 'COR'].copy()
        if not cor_df['atom_name'].equals(itp_df['atomname']):
            sys.exit("The columns name are different!\n")
        cor_charges: pd.DataFrame = itp[['atomnr', 'charge']].copy()
        cor_charges['atomnr'] = cor_charges['atomnr'].astype(str)
        cor_df['atom_id'] = cor_df['atom_id'].astype(str)
        cor_df = pd.merge(cor_df,
                          cor_charges[['atomnr', 'charge']],
                          left_on='atom_id',
                          right_on='atomnr',
                          how='left')
        return cor_df.drop('atomnr', axis=1)

    @staticmethod
    def add_charge_column_to_df(main_df: pd.DataFrame,
                                charge_df: pd.DataFrame
                                ) -> pd.DataFrame:
        """add charge column to the dataframe"""
        tmp_df = main_df.copy()
        tmp_df = pd.merge(tmp_df,
                          charge_df[['atomname', 'charge']],
                          left_on='atom_name',
                          right_on='atomname',
                          how='left')
        return tmp_df.drop('atomname', axis=1)

    @staticmethod
    def get_aptes_charges(itp_df: pd.DataFrame
                          ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """return charges for pro- and unprotonated aptes"""
        pro_found: bool = False
        unpro_found: bool = False
        for _, group in itp_df.groupby('resnr'):
            if len(group) == 13 and not pro_found:
                pro_charges: pd.DataFrame = \
                    group[['atomname', 'charge']].copy()
                pro_found = True
            elif len(group) == 12 and not unpro_found:
                unpro_charges: pd.DataFrame = \
                    group[['atomname', 'charge']].copy()
                unpro_found = True
            if pro_found and unpro_found:
                break
        return pro_charges, unpro_charges


if __name__ == '__main__':
    PdbToPqr(log=logger.setup_logger('pdb2pqr.log'))
