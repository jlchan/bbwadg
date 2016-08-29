function simplex_gm_rule_test01 ( )

%*****************************************************************************80
%
%% SIMPLEX_GM_RULE_TEST01 tests SIMPLEX_UNIT_TO_GENERAL.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 March 2008
%
%  Author:
%
%    John Burkardt
%
  dim_num = 2;
  vertex_num = dim_num + 1;
  point_num = 50;

  t = [
    -1.0, -1.0;
    1.0, -1.0;
    -1.0, 1.0 ]';

  t_unit = [
    0.0, 0.0;
    1.0, 0.0;
    0.0, 1.0 ]';

  seed = 123456789;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'SIMPLEX_GM_RULE_TEST01\n' );
  fprintf ( 1, '  SIMPLEX_UNIT_TO_GENERAL\n' );
  fprintf ( 1, '  maps points in the unit simplex to a general simplex.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Here we consider a simplex in 2D, a triangle.\n' );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  The vertices of the general triangle are:\n' );
  fprintf ( 1, '\n' );
  for j = 1 : vertex_num
    for dim = 1 : dim_num
      fprintf ( 1, '  %8f', t(dim,j) );
    end
    fprintf ( 1, '\n' );
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '   (  XSI     ETA )   ( X       Y  )\n' );
  fprintf ( 1, '\n' );

  phy_unit = simplex_unit_to_general ( dim_num, dim_num + 1, t, t_unit );

  for j = 1 : dim_num + 1
    for dim = 1 : dim_num
      fprintf ( 1, '  %8f', t_unit(dim,j) );
    end
    for dim = 1 : dim_num
      fprintf ( 1, '  %8f', phy_unit(dim,j) );
    end
    fprintf ( 1, '\n' );
  end

  [ ref, seed ] = simplex_unit_sample ( dim_num, point_num, seed );

  phy = simplex_unit_to_general ( dim_num, point_num, t, ref );
% keyboard

plot(phy(1,:),phy(2,:),'o')
  for j = 1 : point_num
    for dim = 1 : dim_num
      fprintf ( 1, '  %8f', ref(dim,j) );
    end
    for dim = 1 : dim_num
      fprintf ( 1, '  %8f', phy(dim,j) );
    end
    fprintf ( 1, '\n' );
  end

  return
end
