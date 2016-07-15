object OptDialog: TOptDialog
  Left = 0
  Top = 0
  BorderIcons = [biSystemMenu]
  BorderStyle = bsDialog
  Caption = 'Options'
  ClientHeight = 644
  ClientWidth = 400
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poDesigned
  OnShow = FormShow
  PixelsPerInch = 96
  TextHeight = 13
  object Label94: TLabel
    Left = 44
    Top = 370
    Width = 149
    Height = 13
    Caption = 'Output Receiver/Satellite Clock'
  end
  object BtnCancel: TButton
    Left = 305
    Top = 607
    Width = 69
    Height = 21
    Caption = '&Cancel'
    ModalResult = 2
    TabOrder = 1
  end
  object BtnOk: TButton
    Left = 235
    Top = 607
    Width = 69
    Height = 21
    Caption = '&OK'
    ModalResult = 1
    TabOrder = 0
    OnClick = BtnOkClick
  end
  object BtnSave: TButton
    Left = 159
    Top = 607
    Width = 69
    Height = 21
    Caption = '&Save'
    TabOrder = 3
    OnClick = BtnSaveClick
  end
  object BtnLoad: TButton
    Left = 84
    Top = 607
    Width = 69
    Height = 21
    Caption = '&Load'
    TabOrder = 2
    OnClick = BtnLoadClick
  end
  object Misc: TPageControl
    Left = 0
    Top = 0
    Width = 400
    Height = 600
    ActivePage = TabSheet7
    Align = alTop
    TabOrder = 4
    object TabSheet1: TTabSheet
      Caption = 'Setting&1'
      object Label3: TLabel
        Left = 20
        Top = 119
        Width = 177
        Height = 13
        Caption = 'Rec Dynamics/Earth Tides Correction'
      end
      object Label8: TLabel
        Left = 20
        Top = 141
        Width = 108
        Height = 13
        Caption = 'Ionosphere Correction'
      end
      object LabelPosMode: TLabel
        Left = 20
        Top = 7
        Width = 80
        Height = 13
        Caption = 'Positioning Mode'
      end
      object LabelFreq: TLabel
        Left = 20
        Top = 29
        Width = 58
        Height = 13
        Caption = 'Frequencies'
      end
      object LabelSolution: TLabel
        Left = 20
        Top = 73
        Width = 65
        Height = 13
        Caption = 'Solution Type'
      end
      object LabelElMask: TLabel
        Left = 20
        Top = 95
        Width = 179
        Height = 13
        Caption = 'Elevation Mask ('#176') / SNR Mask (dbHz)'
      end
      object Label32: TLabel
        Left = 20
        Top = 207
        Width = 119
        Height = 13
        Caption = 'Satellite Ephemeris/Clock'
      end
      object Label35: TLabel
        Left = 20
        Top = 253
        Width = 176
        Height = 13
        Caption = 'Excluded Satellites (+PRN: Included)'
      end
      object Label9: TLabel
        Left = 20
        Top = 163
        Width = 114
        Height = 13
        Caption = 'Troposphere Correction'
      end
      object Label52: TLabel
        Left = 20
        Top = 51
        Width = 76
        Height = 13
        Caption = 'L2 Code Priority'
      end
      object Label53: TLabel
        Left = 20
        Top = 185
        Width = 113
        Height = 13
        Caption = 'Time System Correction'
      end
      object Label87: TLabel
        Left = 20
        Top = 295
        Width = 116
        Height = 13
        Caption = 'Glonass L1 Code Priority'
      end
      object Label88: TLabel
        Left = 20
        Top = 317
        Width = 116
        Height = 13
        Caption = 'Glonass L2 Code Priority'
      end
      object DynamicModel: TComboBox
        Left = 233
        Top = 116
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 2
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'ON')
      end
      object IonoOpt: TComboBox
        Left = 233
        Top = 138
        Width = 113
        Height = 21
        Style = csDropDownList
        DropDownCount = 16
        ItemIndex = 0
        TabOrder = 3
        Text = 'OFF'
        OnChange = IonoOptChange
        Items.Strings = (
          'OFF'
          'Broadcast'
          'SBAS'
          'Dual-Frequency'
          'Estimate STEC'
          'IONEX TEC'
          'QZSS Broardcast')
      end
      object TropOpt: TComboBox
        Left = 233
        Top = 160
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 4
        Text = 'OFF'
        OnChange = TropOptChange
        Items.Strings = (
          'OFF'
          'Saastamoinen'
          'SBAS'
          'Estimate ZTD'
          'Estimate ZTD+Grad')
      end
      object PosMode: TComboBox
        Left = 233
        Top = 3
        Width = 113
        Height = 21
        Style = csDropDownList
        DropDownCount = 10
        ItemIndex = 0
        TabOrder = 0
        Text = 'Single'
        OnChange = PosModeChange
        Items.Strings = (
          'Single'
          'DGPS/DGNSS'
          'Kinematic'
          'Static'
          'Moving-Base'
          'Fixed'
          'PPP Kinematic'
          'PPP Static'
          'PPP Fixed'
          'Multi Baseline Static')
      end
      object Solution: TComboBox
        Left = 233
        Top = 70
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 1
        Text = 'Forward'
        Items.Strings = (
          'Forward'
          'Backward'
          'Combined')
      end
      object SatEphem: TComboBox
        Left = 233
        Top = 204
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 5
        Text = 'Broadcast'
        OnChange = SatEphemChange
        OnClick = SatEphemClick
        Items.Strings = (
          'Broadcast'
          'Precise'
          'Broadcast+SBAS'
          'Broadcast+SSR APC'
          'Broadcast+SSR CoM')
      end
      object ExSats: TEdit
        Left = 233
        Top = 250
        Width = 113
        Height = 21
        TabOrder = 6
      end
      object NavSys1: TCheckBox
        Left = 38
        Top = 275
        Width = 49
        Height = 17
        Caption = 'GPS'
        Checked = True
        State = cbChecked
        TabOrder = 7
      end
      object NavSys2: TCheckBox
        Left = 82
        Top = 275
        Width = 71
        Height = 17
        Caption = 'GLO'
        TabOrder = 8
        OnClick = NavSys2Click
      end
      object NavSys3: TCheckBox
        Left = 128
        Top = 275
        Width = 61
        Height = 17
        Caption = 'Galileo'
        TabOrder = 9
      end
      object NavSys4: TCheckBox
        Left = 182
        Top = 275
        Width = 61
        Height = 17
        Caption = 'QZSS'
        TabOrder = 10
      end
      object NavSys5: TCheckBox
        Left = 234
        Top = 275
        Width = 61
        Height = 17
        Caption = 'SBAS'
        TabOrder = 11
      end
      object TideCorr: TComboBox
        Left = 290
        Top = 116
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 12
        Text = 'OFF'
        OnChange = DynamicModelChange
        Items.Strings = (
          'OFF'
          'ON')
      end
      object NavSys6: TCheckBox
        Left = 286
        Top = 275
        Width = 61
        Height = 17
        Caption = 'Beidou'
        TabOrder = 13
      end
      object ElMask: TComboBox
        Left = 233
        Top = 92
        Width = 56
        Height = 21
        AutoComplete = False
        DropDownCount = 16
        TabOrder = 14
        Text = '15'
        OnChange = DynamicModelChange
        Items.Strings = (
          '0'
          '5'
          '10'
          '15'
          '20'
          '25'
          '30'
          '35'
          '40'
          '45'
          '50'
          '55'
          '60'
          '65'
          '70')
      end
      object L2Cod: TComboBox
        Left = 233
        Top = 48
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 15
        Text = 'L2P(Y)'
        OnChange = FreqChange
        Items.Strings = (
          'L2P(Y)'
          'L2C')
      end
      object TimSys: TComboBox
        Left = 233
        Top = 182
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 16
        Text = 'OFF'
        OnChange = TropOptChange
        Items.Strings = (
          'OFF'
          'Corrections'
          'Nav+EST')
      end
      object BtnMask: TButton
        Left = 290
        Top = 92
        Width = 56
        Height = 22
        Caption = '...'
        TabOrder = 17
        OnClick = BtnMaskClick
      end
      object PosOpt5: TCheckBox
        Left = 298
        Top = 230
        Width = 72
        Height = 17
        Caption = 'RAIM FDE'
        TabOrder = 18
      end
      object PosOpt4: TCheckBox
        Left = 226
        Top = 231
        Width = 66
        Height = 17
        Caption = 'Reject Ecl'
        TabOrder = 19
      end
      object PosOpt3: TCheckBox
        Left = 157
        Top = 230
        Width = 69
        Height = 17
        Caption = 'PhWindup'
        TabOrder = 20
      end
      object PosOpt2: TCheckBox
        Left = 80
        Top = 230
        Width = 57
        Height = 17
        Caption = 'Rec PCV'
        TabOrder = 21
      end
      object PosOpt1: TCheckBox
        Left = 3
        Top = 230
        Width = 62
        Height = 17
        Caption = 'Sat PCV'
        TabOrder = 22
      end
      object Freqs: TComboBox
        Left = 233
        Top = 26
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 23
        Text = 'L1'
        OnChange = FreqChange
        Items.Strings = (
          'L1'
          'L5'
          'L1+L2'
          'L1+L5'
          'L1+L2+L5')
      end
      object GloCodePri1: TEdit
        Left = 233
        Top = 294
        Width = 113
        Height = 21
        TabOrder = 24
        Text = 'PC'
        OnKeyPress = GloCodePri1KeyPress
      end
      object GloCodePri2: TEdit
        Left = 233
        Top = 316
        Width = 113
        Height = 21
        TabOrder = 25
        Text = 'PC'
        OnKeyPress = GloCodePri2KeyPress
      end
    end
    object TabSheet2: TTabSheet
      Caption = 'Setting&2'
      ImageIndex = 1
      object Label25: TLabel
        Left = 20
        Top = 3
        Width = 178
        Height = 13
        Caption = 'Integer Ambiguity Resolution Method'
        Layout = tlCenter
      end
      object Label24: TLabel
        Left = 20
        Top = 91
        Width = 124
        Height = 13
        Caption = 'Min Ratio to Fix Ambiguity'
        Layout = tlCenter
      end
      object Label13: TLabel
        Left = 20
        Top = 135
        Width = 190
        Height = 13
        Caption = 'Min Lock / Elevation ('#176') to Fix Ambiguity'
        Layout = tlCenter
      end
      object LabelHold: TLabel
        Left = 20
        Top = 157
        Width = 190
        Height = 13
        Caption = 'Min Fix / Elevation ('#176') to Hold Ambiguity'
        Layout = tlCenter
      end
      object Label22: TLabel
        Left = 20
        Top = 179
        Width = 173
        Height = 13
        Caption = 'Outage to Reset Amb/Slip Thres (m)'
        Layout = tlCenter
      end
      object Label14: TLabel
        Left = 20
        Top = 245
        Width = 127
        Height = 13
        Caption = 'Max Age of Differential (s)'
        Layout = tlCenter
      end
      object Label11: TLabel
        Left = 20
        Top = 267
        Width = 176
        Height = 13
        Caption = 'Reject Threshold of GDOP/Innov (m)'
        Layout = tlCenter
      end
      object Label37: TLabel
        Left = 20
        Top = 289
        Width = 122
        Height = 13
        Caption = 'Number of Filter Iteration'
        Layout = tlCenter
      end
      object Label33: TLabel
        Left = 20
        Top = 47
        Width = 149
        Height = 13
        Caption = 'GLONASS Ambiguity Resolution'
        Layout = tlCenter
      end
      object Label55: TLabel
        Left = 20
        Top = 201
        Width = 83
        Height = 13
        Caption = 'Phase Cycle Shift'
        Layout = tlCenter
      end
      object Label67: TLabel
        Left = 20
        Top = 333
        Width = 84
        Height = 13
        Caption = 'Inter System Bias'
        Layout = tlCenter
      end
      object Label69: TLabel
        Left = 20
        Top = 355
        Width = 190
        Height = 13
        Caption = 'Analysys Method in Double Differencing'
        Layout = tlCenter
      end
      object Label77: TLabel
        Left = 20
        Top = 69
        Width = 121
        Height = 13
        Caption = 'PPP Ambiguity Resolution'
        Layout = tlCenter
      end
      object Label80: TLabel
        Left = 20
        Top = 223
        Width = 61
        Height = 13
        Caption = 'L2C-L2P Bias'
        Layout = tlCenter
      end
      object LabelConf: TLabel
        Left = 20
        Top = 113
        Width = 179
        Height = 13
        Caption = 'Min Confidence / Max FCB to Fix Amb'
        Layout = tlCenter
      end
      object Label84: TLabel
        Left = 20
        Top = 25
        Width = 184
        Height = 13
        Caption = 'Integer Ambiguity Resolution Strategy'
        Layout = tlCenter
      end
      object AmbRes: TComboBox
        Left = 233
        Top = 25
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 0
        Text = 'Continuous'
        OnChange = AmbResChange
        Items.Strings = (
          'Continuous'
          'Instantaneous'
          'Fix and Hold')
      end
      object ValidThresAR: TEdit
        Left = 233
        Top = 91
        Width = 113
        Height = 21
        TabOrder = 2
        Text = '3.0'
      end
      object LockCntFixAmb: TEdit
        Left = 233
        Top = 135
        Width = 56
        Height = 21
        TabOrder = 4
        Text = '5'
      end
      object OutCntResetAmb: TEdit
        Left = 233
        Top = 179
        Width = 56
        Height = 21
        TabOrder = 6
        Text = '5'
      end
      object ElMaskAR: TEdit
        Left = 290
        Top = 135
        Width = 56
        Height = 21
        TabOrder = 3
        Text = '0'
      end
      object SlipThres: TEdit
        Left = 290
        Top = 179
        Width = 56
        Height = 21
        TabOrder = 7
        Text = '0.05'
      end
      object MaxAgeDiff: TEdit
        Left = 233
        Top = 245
        Width = 113
        Height = 21
        TabOrder = 8
        Text = '30'
      end
      object RejectThres: TEdit
        Left = 290
        Top = 267
        Width = 56
        Height = 21
        TabOrder = 10
        Text = '30'
      end
      object NumIter: TEdit
        Left = 233
        Top = 289
        Width = 113
        Height = 21
        TabOrder = 11
        Text = '1'
      end
      object BaselineLen: TEdit
        Left = 233
        Top = 311
        Width = 56
        Height = 21
        TabOrder = 13
        Text = '0.0'
      end
      object BaselineSig: TEdit
        Left = 290
        Top = 311
        Width = 56
        Height = 21
        TabOrder = 14
        Text = '0.001'
      end
      object BaselineConst: TCheckBox
        Left = 20
        Top = 311
        Width = 179
        Height = 21
        Caption = 'Baseline Length Constraint (m)'
        TabOrder = 12
        OnClick = BaselineConstClick
      end
      object GloAmbRes: TComboBox
        Left = 233
        Top = 47
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 1
        Text = 'OFF'
        OnChange = AmbResChange
        Items.Strings = (
          'OFF'
          'ON'
          'Auto Calibration'
          'Use IFB Table')
      end
      object FixCntHoldAmb: TEdit
        Left = 233
        Top = 157
        Width = 56
        Height = 21
        TabOrder = 15
        Text = '10'
      end
      object ElMaskHold: TEdit
        Left = 290
        Top = 157
        Width = 56
        Height = 21
        TabOrder = 5
        Text = '0'
      end
      object RejectGdop: TEdit
        Left = 233
        Top = 267
        Width = 56
        Height = 21
        TabOrder = 9
        Text = '30'
      end
      object PhaCyc: TComboBox
        Left = 233
        Top = 201
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 16
        Text = 'OFF'
        OnChange = AmbResChange
        Items.Strings = (
          'OFF'
          'Table'
          'RINEX/RTCM')
      end
      object Isb: TComboBox
        Left = 233
        Top = 333
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 17
        Text = 'OFF'
        OnChange = IsbChange
        Items.Strings = (
          'OFF'
          'Table'
          'Estimate Phase+Code'
          'Estimate Code'
          'Estimate Phase'
          'Estimation(0m BL)')
      end
      object Diff: TComboBox
        Left = 233
        Top = 355
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 18
        Text = 'inall'
        OnChange = AmbResChange
        Items.Strings = (
          'inall'
          'exc. glonass')
      end
      object PppAmbRes: TComboBox
        Left = 233
        Top = 69
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 19
        Text = 'OFF'
        OnChange = AmbResChange
        Items.Strings = (
          'OFF'
          'CNES'
          'CNES-ILS'
          'FCB')
      end
      object L2CPBias: TComboBox
        Left = 233
        Top = 223
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 20
        Text = 'OFF'
        OnChange = L2CPBiasChange
        Items.Strings = (
          'OFF'
          'Table'
          'Estimation')
      end
      object ThresAR2: TEdit
        Left = 233
        Top = 113
        Width = 56
        Height = 21
        TabOrder = 21
        Text = '0.99995'
      end
      object ThresAR3: TEdit
        Left = 290
        Top = 113
        Width = 56
        Height = 21
        TabOrder = 22
        Text = '0.20'
      end
      object AmbResMethod: TComboBox
        Left = 233
        Top = 3
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 23
        Text = 'OFF'
        OnChange = AmbResChange
        Items.Strings = (
          'OFF'
          'LAMBDA'
          'WL/NL'
          'TCAR')
      end
    end
    object TabSheet8: TTabSheet
      Caption = 'Setting&3'
      ImageIndex = 2
      object Label12: TLabel
        Left = 4
        Top = 7
        Width = 217
        Height = 13
        Caption = 'Phase Cycle Shift, GLONASS IFB, Error Model'
      end
      object PhaCycSBtn: TSpeedButton
        Left = 313
        Top = 5
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = PhaCycSBtnClick
      end
      object GloIfbSBtn: TSpeedButton
        Left = 333
        Top = 5
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = GloIfbSBtnClick
      end
      object ErrModSBtn: TSpeedButton
        Left = 354
        Top = 5
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = ErrModSBtnClick
      end
      object PhaCycBtn: TButton
        Left = 354
        Top = 25
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 0
        OnClick = PhaCycBtnClick
      end
      object PhaCycFile: TEdit
        Left = 0
        Top = 22
        Width = 353
        Height = 21
        TabOrder = 1
      end
      object GloIfbFile: TEdit
        Left = 0
        Top = 44
        Width = 353
        Height = 21
        TabOrder = 2
      end
      object ErrModFile: TEdit
        Left = 0
        Top = 66
        Width = 353
        Height = 21
        TabOrder = 3
      end
      object ErrModBtn: TButton
        Left = 354
        Top = 69
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 4
        OnClick = ErrModBtnClick
      end
      object GloIfbBtn: TButton
        Left = 354
        Top = 47
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 5
        OnClick = GloIfbBtnClick
      end
      object GroupBox1: TGroupBox
        Left = -1
        Top = 92
        Width = 371
        Height = 389
        Caption = 'Multi Baseline Static'
        TabOrder = 6
        object Label49: TLabel
          Left = 13
          Top = 20
          Width = 133
          Height = 13
          Caption = 'Estimate Satellite Clock/FCB'
        end
        object Label50: TLabel
          Left = 13
          Top = 41
          Width = 172
          Height = 13
          Caption = 'Semi-Dynamic Correction Parameter'
        end
        object Label17: TLabel
          Left = 13
          Top = 79
          Width = 85
          Height = 13
          Caption = 'Solution Directory'
        end
        object Label68: TLabel
          Left = 13
          Top = 206
          Width = 115
          Height = 13
          Caption = 'Fixing Probability WL/NL'
        end
        object Label66: TLabel
          Left = 13
          Top = 188
          Width = 152
          Height = 13
          Caption = 'O-C Reject Phase/Code (sigma)'
        end
        object Label61: TLabel
          Left = 13
          Top = 167
          Width = 153
          Height = 13
          Caption = 'Trop. Process Noise Zen/EW/NS'
        end
        object Label60: TLabel
          Left = 13
          Top = 146
          Width = 162
          Height = 13
          Caption = 'Est. Interval of Trop. Gradient (s)'
        end
        object Label48: TLabel
          Left = 13
          Top = 125
          Width = 111
          Height = 13
          Caption = 'Est. Interval of ZTD (s)'
        end
        object Label70: TLabel
          Left = 13
          Top = 224
          Width = 156
          Height = 13
          Caption = 'Convergence Factor of Iteration'
        end
        object Label71: TLabel
          Left = 13
          Top = 245
          Width = 89
          Height = 13
          Caption = 'Maximum Iteration'
        end
        object Label54: TLabel
          Left = 13
          Top = 263
          Width = 188
          Height = 13
          Caption = 'Temporal Storage of Epoch Parameters'
        end
        object Label23: TLabel
          Left = 13
          Top = 306
          Width = 151
          Height = 13
          Caption = 'Minimum Pass Length ZD/DD (s)'
        end
        object Label78: TLabel
          Left = 13
          Top = 325
          Width = 157
          Height = 13
          Caption = 'Mobile Station STD dN/dE/dU (m)'
        end
        object Label79: TLabel
          Left = 13
          Top = 347
          Width = 150
          Height = 13
          Caption = 'Base Station STD dN/dE/dU (m)'
        end
        object Label93: TLabel
          Left = 13
          Top = 370
          Width = 123
          Height = 13
          Caption = 'NL FCB time interval (sec)'
        end
        object EstSatClo: TComboBox
          Left = 229
          Top = 17
          Width = 56
          Height = 21
          Style = csDropDownList
          ItemIndex = 0
          TabOrder = 0
          Text = 'OFF'
          Items.Strings = (
            'OFF'
            'ON')
        end
        object EstSatFCB: TComboBox
          Left = 291
          Top = 17
          Width = 56
          Height = 21
          Style = csDropDownList
          ItemIndex = 0
          TabOrder = 1
          Text = 'OFF'
          OnChange = EstSatFCBChange
          Items.Strings = (
            'OFF'
            'ON')
        end
        object SemiDCPara: TEdit
          Left = 13
          Top = 55
          Width = 334
          Height = 21
          TabOrder = 2
        end
        object SemiDCParaBtn: TButton
          Left = 350
          Top = 58
          Width = 17
          Height = 17
          Caption = '...'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -9
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 3
          OnClick = SemiDCParaBtnClick
        end
        object SolDirFile: TEdit
          Left = 13
          Top = 95
          Width = 334
          Height = 21
          TabOrder = 4
        end
        object SolDirBtn: TButton
          Left = 351
          Top = 99
          Width = 17
          Height = 17
          Caption = '...'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -9
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 5
          OnClick = SolDirBtnClick
        end
        object EstIntES: TEdit
          Left = 226
          Top = 143
          Width = 121
          Height = 21
          TabOrder = 6
          Text = '43200'
        end
        object RanWalSigZen: TEdit
          Left = 226
          Top = 164
          Width = 40
          Height = 21
          TabOrder = 7
          Text = '1.00E-04'
        end
        object ThrO: TEdit
          Left = 226
          Top = 185
          Width = 60
          Height = 21
          TabOrder = 8
          Text = '5.0'
        end
        object EstIntZ: TEdit
          Left = 226
          Top = 122
          Width = 121
          Height = 21
          TabOrder = 9
          Text = '7200'
        end
        object ConCri: TEdit
          Left = 226
          Top = 221
          Width = 121
          Height = 21
          TabOrder = 10
          Text = '1.00E-03'
        end
        object MaxIte: TEdit
          Left = 226
          Top = 242
          Width = 121
          Height = 21
          TabOrder = 11
          Text = '3'
        end
        object RanWalSigEW: TEdit
          Left = 266
          Top = 164
          Width = 41
          Height = 21
          TabOrder = 12
          Text = '1.00E-04'
        end
        object RanWalSigNS: TEdit
          Left = 307
          Top = 164
          Width = 40
          Height = 21
          TabOrder = 13
          Text = '1.00E-07'
        end
        object ThrC: TEdit
          Left = 286
          Top = 185
          Width = 61
          Height = 21
          TabOrder = 14
          Text = '5.0'
        end
        object JudValWL: TEdit
          Left = 226
          Top = 203
          Width = 60
          Height = 21
          TabOrder = 15
          Text = '0.99990'
        end
        object JudValL1: TEdit
          Left = 286
          Top = 203
          Width = 61
          Height = 21
          TabOrder = 16
          Text = '0.99990'
        end
        object TemStoFile: TEdit
          Left = 13
          Top = 279
          Width = 334
          Height = 21
          TabOrder = 17
        end
        object TemStoFileBut: TButton
          Left = 350
          Top = 282
          Width = 17
          Height = 17
          Caption = '...'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -9
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 18
          OnClick = TemStoFileButClick
        end
        object Minsd: TEdit
          Left = 226
          Top = 303
          Width = 60
          Height = 21
          TabOrder = 19
          Text = '600'
        end
        object Mindd: TEdit
          Left = 287
          Top = 303
          Width = 60
          Height = 21
          TabOrder = 20
          Text = '1200'
        end
        object Mobstadn: TEdit
          Left = 226
          Top = 322
          Width = 40
          Height = 21
          TabOrder = 21
          Text = '0'
        end
        object Mobstade: TEdit
          Left = 267
          Top = 322
          Width = 40
          Height = 21
          TabOrder = 22
          Text = '0'
        end
        object Mobstadu: TEdit
          Left = 307
          Top = 322
          Width = 40
          Height = 21
          TabOrder = 23
          Text = '0'
        end
        object Basestadu: TEdit
          Left = 307
          Top = 344
          Width = 40
          Height = 21
          TabOrder = 24
          Text = '0'
        end
        object Basestade: TEdit
          Left = 267
          Top = 344
          Width = 40
          Height = 21
          TabOrder = 25
          Text = '0'
        end
        object Basestadn: TEdit
          Left = 226
          Top = 344
          Width = 40
          Height = 21
          TabOrder = 26
          Text = '0'
        end
        object NLFcb: TEdit
          Left = 226
          Top = 368
          Width = 121
          Height = 21
          TabOrder = 27
          Text = '7200'
        end
      end
      object GroupBox2: TGroupBox
        Left = 0
        Top = 487
        Width = 371
        Height = 90
        Caption = 'PPP-AR'
        TabOrder = 7
        object FCBLabel: TLabel
          Left = 13
          Top = 18
          Width = 60
          Height = 13
          Caption = 'Satellite FCB'
          Enabled = False
        end
        object Label75: TLabel
          Left = 14
          Top = 62
          Width = 115
          Height = 13
          Caption = 'Fixing Probability WL/NL'
        end
        object FCBFile: TEdit
          Left = 13
          Top = 37
          Width = 334
          Height = 21
          TabOrder = 0
        end
        object FCBBtn: TButton
          Left = 350
          Top = 40
          Width = 17
          Height = 17
          Caption = '...'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -9
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 1
          OnClick = FCBBtnClick
        end
        object Fixnl: TEdit
          Left = 288
          Top = 61
          Width = 60
          Height = 21
          TabOrder = 2
          Text = '0.99990'
        end
        object Fixwl: TEdit
          Left = 227
          Top = 61
          Width = 60
          Height = 21
          TabOrder = 3
          Text = '0.99990'
        end
      end
    end
    object TabSheet3: TTabSheet
      Caption = 'O&utput'
      ImageIndex = 3
      object LabelSolFormat: TLabel
        Left = 36
        Top = 6
        Width = 75
        Height = 13
        Caption = 'Solution Format'
        Layout = tlCenter
      end
      object LabelTimeFormat: TLabel
        Left = 36
        Top = 50
        Width = 134
        Height = 13
        Caption = 'Time Format / # of Decimals'
        Layout = tlCenter
      end
      object LabelLatLonFormat: TLabel
        Left = 36
        Top = 72
        Width = 133
        Height = 13
        Caption = 'Latitude / Longitude Format'
        Layout = tlCenter
      end
      object LabelFieldSep: TLabel
        Left = 36
        Top = 94
        Width = 73
        Height = 13
        Caption = 'Field Separator'
        Layout = tlCenter
      end
      object Label2: TLabel
        Left = 36
        Top = 116
        Width = 66
        Height = 13
        Caption = 'Datum/Height'
        Layout = tlCenter
      end
      object Label18: TLabel
        Left = 36
        Top = 138
        Width = 58
        Height = 13
        Caption = 'Geoid Model'
        Layout = tlCenter
      end
      object Label20: TLabel
        Left = 36
        Top = 28
        Width = 167
        Height = 13
        Caption = 'Output Header/Processing Options'
        Layout = tlCenter
      end
      object Label36: TLabel
        Left = 36
        Top = 204
        Width = 180
        Height = 13
        Caption = 'Output Solution Status / Debug Trace'
        Layout = tlCenter
      end
      object Label21: TLabel
        Left = 36
        Top = 182
        Width = 185
        Height = 13
        Caption = 'NMEA Interval (s) RMC/GGA, GSA/GSV'
        Enabled = False
        Layout = tlCenter
      end
      object Label31: TLabel
        Left = 36
        Top = 160
        Width = 114
        Height = 13
        Caption = 'Solution for Static Mode'
        Layout = tlCenter
      end
      object Label76: TLabel
        Left = 36
        Top = 226
        Width = 79
        Height = 13
        Caption = 'Output ISB Data'
        Layout = tlCenter
      end
      object BtnISBOutView: TSpeedButton
        Left = 347
        Top = 228
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnISBOutViewClick
      end
      object Label89: TLabel
        Left = 36
        Top = 362
        Width = 149
        Height = 13
        Caption = 'Output Receiver/Satellite Clock'
      end
      object Label90: TLabel
        Left = 36
        Top = 340
        Width = 79
        Height = 13
        Caption = 'Output Ion/Trop'
      end
      object Label91: TLabel
        Left = 36
        Top = 318
        Width = 117
        Height = 13
        Caption = 'Output Position in SINEX'
      end
      object Label92: TLabel
        Left = 36
        Top = 274
        Width = 102
        Height = 13
        Caption = 'Output L2P-L2C Data'
        Layout = tlCenter
      end
      object BtnGL2OutView: TSpeedButton
        Left = 347
        Top = 276
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnGL2OutViewClick
      end
      object Label95: TLabel
        Left = 36
        Top = 384
        Width = 76
        Height = 13
        Caption = 'Output Baseline'
      end
      object Label96: TLabel
        Left = 36
        Top = 406
        Width = 56
        Height = 13
        Caption = 'Output FCB'
      end
      object BtnFCBOutView: TSpeedButton
        Left = 347
        Top = 407
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnFCBOutViewClick
      end
      object SolFormat: TComboBox
        Left = 233
        Top = 3
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 0
        Text = 'Lat/Lon/Height'
        OnChange = SolFormatChange
        Items.Strings = (
          'Lat/Lon/Height'
          'X/Y/Z-ECEF'
          'E/N/U-Baseline'
          'NMEA0183'
          'R/V/A-ECEF')
      end
      object TimeFormat: TComboBox
        Left = 233
        Top = 47
        Width = 96
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 3
        Text = 'ww ssss GPST'
        Items.Strings = (
          'ww ssss GPST'
          'hh:mm:ss GPST'
          'hh:mm:ss UTC'
          'hh:mm:ss JST')
      end
      object LatLonFormat: TComboBox
        Left = 233
        Top = 69
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 4
        Text = 'ddd.ddddddd'
        Items.Strings = (
          'ddd.ddddddd'
          'ddd mm ss.sss')
      end
      object FieldSep: TEdit
        Left = 233
        Top = 91
        Width = 113
        Height = 21
        TabOrder = 5
      end
      object OutputDatum: TComboBox
        Left = 233
        Top = 113
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 6
        Text = 'WGS84'
        Items.Strings = (
          'WGS84')
      end
      object OutputHeight: TComboBox
        Left = 290
        Top = 113
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 7
        Text = 'Ellipsoidal'
        OnClick = OutputHeightClick
        Items.Strings = (
          'Ellipsoidal'
          'Geodetic')
      end
      object OutputGeoid: TComboBox
        Left = 233
        Top = 135
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 8
        Text = 'Internal'
        Items.Strings = (
          'Internal'
          'EGM96-BE (15")'
          'EGM2008-SE (2.5")'
          'EGM2008-SE (1.0")'
          'GSI2000 (1x1.5")')
      end
      object OutputHead: TComboBox
        Left = 233
        Top = 25
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 1
        TabOrder = 1
        Text = 'ON'
        Items.Strings = (
          'OFF'
          'ON')
      end
      object OutputOpt: TComboBox
        Left = 290
        Top = 25
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 2
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'ON')
      end
      object DebugTrace: TComboBox
        Left = 290
        Top = 201
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 9
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'Level1'
          'Level2'
          'Level3'
          'Level4'
          'Level5')
      end
      object DebugStatus: TComboBox
        Left = 233
        Top = 201
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 10
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'State'
          'Residuals')
      end
      object TimeDecimal: TEdit
        Left = 329
        Top = 47
        Width = 17
        Height = 21
        TabOrder = 11
        Text = '3'
      end
      object NmeaIntv1: TEdit
        Left = 233
        Top = 179
        Width = 56
        Height = 21
        Enabled = False
        TabOrder = 12
        Text = '0'
      end
      object NmeaIntv2: TEdit
        Left = 290
        Top = 179
        Width = 56
        Height = 21
        Enabled = False
        TabOrder = 13
        Text = '0'
      end
      object SolStatic: TComboBox
        Left = 233
        Top = 157
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 14
        Text = 'All'
        Items.Strings = (
          'All'
          'Single')
      end
      object ISBOutFile: TEdit
        Left = 46
        Top = 245
        Width = 300
        Height = 21
        TabOrder = 15
      end
      object BtnISBOutFile: TButton
        Left = 347
        Top = 247
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 16
        OnClick = BtnISBOutFileClick
      end
      object IsbOut: TComboBox
        Left = 233
        Top = 223
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 17
        Text = 'OFF'
        OnChange = IsbOutChange
        Items.Strings = (
          'OFF'
          'NEW'
          'APPEND')
      end
      object OutSatelliteClock: TComboBox
        Left = 290
        Top = 355
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 18
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'OUT')
      end
      object OutReceiverClock: TComboBox
        Left = 233
        Top = 355
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 19
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'OUT')
      end
      object OutPosSinex: TComboBox
        Left = 233
        Top = 311
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 20
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'OUT')
      end
      object OutIon: TComboBox
        Left = 233
        Top = 333
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 21
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'OUT')
      end
      object OutTrop: TComboBox
        Left = 290
        Top = 333
        Width = 56
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 22
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'OUT')
      end
      object Gl2Out: TComboBox
        Left = 233
        Top = 267
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 23
        Text = 'OFF'
        OnChange = Gl2OutChange
        Items.Strings = (
          'OFF'
          'NEW'
          'APPEND')
      end
      object BtnGL2OutFile: TButton
        Left = 347
        Top = 293
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 24
        OnClick = BtnGL2OutFileClick
      end
      object GL2OutFile: TEdit
        Left = 46
        Top = 289
        Width = 300
        Height = 21
        TabOrder = 25
      end
      object OutStatic: TComboBox
        Left = 233
        Top = 377
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 26
        Text = 'OFF'
        Items.Strings = (
          'OFF'
          'OUT')
      end
      object FCBOutFile: TEdit
        Left = 46
        Top = 421
        Width = 300
        Height = 21
        TabOrder = 27
      end
      object BtnFCBOutFile: TButton
        Left = 347
        Top = 424
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 28
        OnClick = BtnFCBOutFileClick
      end
      object FCBOut: TComboBox
        Left = 233
        Top = 399
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 29
        Text = 'OFF'
        OnChange = FCBOutChange
        Items.Strings = (
          'OFF'
          'NEW')
      end
    end
    object TabSheet4: TTabSheet
      Caption = 'S&tatistics'
      ImageIndex = 4
      object Label29: TLabel
        Left = 36
        Top = 403
        Width = 132
        Height = 13
        Caption = 'Satellite Clock Stability (s/s)'
        Layout = tlCenter
      end
      object GroupBox3: TGroupBox
        Left = 2
        Top = 0
        Width = 369
        Height = 201
        Caption = 'Measurement Errors (1-sigma)'
        TabOrder = 0
        object Label6: TLabel
          Left = 34
          Top = 42
          Width = 164
          Height = 13
          Caption = 'Code/Carrier-Phase Error Ratio L1'
          Layout = tlCenter
        end
        object Label7: TLabel
          Left = 34
          Top = 108
          Width = 160
          Height = 13
          Caption = 'Carrier-Phase Error a+b/sinEl (m)'
          Layout = tlCenter
        end
        object Label16: TLabel
          Left = 34
          Top = 152
          Width = 184
          Height = 13
          Caption = 'Carrier-Phase Error/Baseline (m/10km)'
          Layout = tlCenter
        end
        object Label64: TLabel
          Left = 34
          Top = 174
          Width = 114
          Height = 13
          Caption = 'Doppler Frequency (Hz)'
          Layout = tlCenter
        end
        object Label56: TLabel
          Left = 34
          Top = 23
          Width = 55
          Height = 13
          Caption = 'Error Model'
          Layout = tlCenter
        end
        object Label57: TLabel
          Left = 34
          Top = 130
          Width = 126
          Height = 13
          Caption = 'Code Error Ratio (no DCB)'
          Layout = tlCenter
        end
        object Label34: TLabel
          Left = 34
          Top = 64
          Width = 164
          Height = 13
          Caption = 'Code/Carrier-Phase Error Ratio L2'
          Layout = tlCenter
        end
        object Label83: TLabel
          Left = 34
          Top = 86
          Width = 164
          Height = 13
          Caption = 'Code/Carrier-Phase Error Ratio L5'
          Layout = tlCenter
        end
        object MeasErrR1: TEdit
          Left = 231
          Top = 42
          Width = 113
          Height = 21
          TabOrder = 0
          Text = '100.0'
        end
        object MeasErr2: TEdit
          Left = 231
          Top = 108
          Width = 56
          Height = 21
          TabOrder = 1
          Text = '0.003'
        end
        object MeasErr3: TEdit
          Left = 288
          Top = 108
          Width = 56
          Height = 21
          TabOrder = 2
          Text = '0.003'
        end
        object MeasErr4: TEdit
          Left = 231
          Top = 152
          Width = 113
          Height = 21
          TabOrder = 3
          Text = '0.000'
        end
        object MeasErr5: TEdit
          Left = 231
          Top = 174
          Width = 113
          Height = 21
          TabOrder = 4
          Text = '0.100'
        end
        object MeasErrR2: TEdit
          Left = 231
          Top = 64
          Width = 113
          Height = 21
          TabOrder = 5
          Text = '100.0'
        end
        object ErrMod: TComboBox
          Left = 231
          Top = 20
          Width = 113
          Height = 21
          Style = csDropDownList
          ItemIndex = 0
          TabOrder = 6
          Text = 'User Settings'
          OnChange = AmbResChange
          Items.Strings = (
            'User Settings'
            'Table')
        end
        object MeasErr6: TEdit
          Left = 231
          Top = 130
          Width = 113
          Height = 21
          TabOrder = 7
          Text = '10.0'
        end
        object MeasErrR3: TEdit
          Left = 231
          Top = 86
          Width = 113
          Height = 21
          TabOrder = 8
          Text = '100.0'
        end
      end
      object GroupBox4: TGroupBox
        Left = 2
        Top = 215
        Width = 369
        Height = 179
        Caption = 'Process Noises (1-sigma/sqrt(s))'
        TabOrder = 1
        object Label26: TLabel
          Left = 34
          Top = 42
          Width = 123
          Height = 13
          Caption = 'Carrier-Phase Bias (cycle)'
          Layout = tlCenter
        end
        object Label27: TLabel
          Left = 34
          Top = 64
          Width = 172
          Height = 13
          Caption = 'Vertical Ionospheric Delay (m/10km)'
          Layout = tlCenter
        end
        object Label28: TLabel
          Left = 34
          Top = 86
          Width = 144
          Height = 13
          Caption = 'Zenith Tropospheric Delay (m)'
          Layout = tlCenter
        end
        object Label10: TLabel
          Left = 34
          Top = 22
          Width = 170
          Height = 13
          Caption = 'Receiver Accel Horiz/Vertical (m/s2)'
          Layout = tlCenter
        end
        object Label72: TLabel
          Left = 34
          Top = 108
          Width = 173
          Height = 13
          Caption = 'Carrier-Phase Inter-System Bias (m)'
          Layout = tlCenter
        end
        object Label73: TLabel
          Left = 34
          Top = 130
          Width = 170
          Height = 13
          Caption = 'Pseudorange Inter-System Bias (m)'
          Layout = tlCenter
        end
        object Label74: TLabel
          Left = 34
          Top = 152
          Width = 80
          Height = 13
          Caption = 'L2C-L2P Bias (m)'
          Layout = tlCenter
        end
        object PrNoise1: TEdit
          Left = 231
          Top = 42
          Width = 113
          Height = 21
          TabOrder = 2
          Text = '1.0E-04'
        end
        object PrNoise2: TEdit
          Left = 231
          Top = 64
          Width = 113
          Height = 21
          TabOrder = 3
          Text = '1.0E-03'
        end
        object PrNoise3: TEdit
          Left = 231
          Top = 86
          Width = 113
          Height = 21
          TabOrder = 4
          Text = '1.0E-04'
        end
        object PrNoise4: TEdit
          Left = 231
          Top = 20
          Width = 56
          Height = 21
          TabOrder = 0
          Text = '1.0E-04'
        end
        object PrNoise5: TEdit
          Left = 288
          Top = 20
          Width = 56
          Height = 21
          TabOrder = 1
          Text = '1.0E-04'
        end
        object PrNoise6: TEdit
          Left = 231
          Top = 108
          Width = 113
          Height = 21
          TabOrder = 5
          Text = '1.0E-04'
        end
        object PrNoise7: TEdit
          Left = 231
          Top = 130
          Width = 113
          Height = 21
          TabOrder = 6
          Text = '1.0E-04'
        end
        object PrNoise8: TEdit
          Left = 231
          Top = 152
          Width = 113
          Height = 21
          TabOrder = 7
          Text = '1.0E-04'
        end
      end
      object SatClkStab: TEdit
        Left = 233
        Top = 403
        Width = 113
        Height = 21
        TabOrder = 2
        Text = '5.0E-12'
      end
    end
    object TabSheet5: TTabSheet
      Caption = '&Positions'
      ImageIndex = 5
      object Label4: TLabel
        Left = 12
        Top = 12
        Width = 3
        Height = 13
      end
      object Label30: TLabel
        Left = 10
        Top = 252
        Width = 93
        Height = 13
        Caption = 'Station Position File'
      end
      object BtnStaPosView: TSpeedButton
        Left = 336
        Top = 269
        Width = 17
        Height = 17
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnStaPosViewClick
      end
      object GroupRovAnt: TGroupBox
        Left = 3
        Top = 0
        Width = 369
        Height = 126
        Caption = 'Rover'
        TabOrder = 0
        object LabelRovAntD: TLabel
          Left = 210
          Top = 58
          Width = 76
          Height = 13
          Caption = 'Delta-E/N/U (m)'
        end
        object Label58: TLabel
          Left = 8
          Top = 101
          Width = 69
          Height = 13
          Caption = 'Receiver Type'
        end
        object RovAntE: TEdit
          Left = 208
          Top = 74
          Width = 51
          Height = 21
          TabOrder = 7
          Text = '0'
        end
        object RovAntN: TEdit
          Left = 260
          Top = 74
          Width = 51
          Height = 21
          TabOrder = 8
          Text = '0'
        end
        object RovAntU: TEdit
          Left = 312
          Top = 74
          Width = 51
          Height = 21
          TabOrder = 9
          Text = '0'
        end
        object RovPos1: TEdit
          Left = 6
          Top = 36
          Width = 119
          Height = 21
          TabOrder = 1
          Text = '0'
        end
        object RovPos2: TEdit
          Left = 124
          Top = 36
          Width = 119
          Height = 21
          TabOrder = 2
          Text = '0'
        end
        object RovPos3: TEdit
          Left = 244
          Top = 36
          Width = 119
          Height = 21
          TabOrder = 3
          Text = '0'
        end
        object BtnRovPos: TButton
          Left = 345
          Top = 16
          Width = 17
          Height = 17
          Caption = '...'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -9
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 4
          OnClick = BtnRovPosClick
        end
        object RovAntPcv: TCheckBox
          Left = 6
          Top = 58
          Width = 161
          Height = 17
          Caption = 'Antenna Type (*: Auto)'
          TabOrder = 5
          OnClick = RovAntPcvClick
        end
        object RovAnt: TComboBox
          Left = 4
          Top = 74
          Width = 203
          Height = 21
          DropDownCount = 16
          TabOrder = 6
          OnClick = RovAntClick
        end
        object RovPosType: TComboBox
          Left = 6
          Top = 14
          Width = 137
          Height = 21
          Style = csDropDownList
          Enabled = False
          ItemIndex = 0
          TabOrder = 0
          Text = 'Lat/Lon/Height (deg/m)'
          OnChange = RovPosTypeChange
          Items.Strings = (
            'Lat/Lon/Height (deg/m)'
            'Lat/Lon/Height (dms/m)'
            'X/Y/Z-ECEF (m)'
            'Average of Single Pos'
            'Get from Position File'
            'RINEX Header Position')
        end
        object RovRecTyp: TEdit
          Left = 124
          Top = 98
          Width = 238
          Height = 21
          TabOrder = 10
        end
      end
      object GroupRefAnt: TGroupBox
        Left = 3
        Top = 125
        Width = 369
        Height = 125
        Caption = 'Base Station'
        TabOrder = 1
        object LabelRefAntD: TLabel
          Left = 210
          Top = 58
          Width = 76
          Height = 13
          Caption = 'Delta-E/N/U (m)'
        end
        object Label59: TLabel
          Left = 8
          Top = 101
          Width = 69
          Height = 13
          Caption = 'Receiver Type'
        end
        object RefAntE: TEdit
          Left = 208
          Top = 74
          Width = 51
          Height = 21
          TabOrder = 7
          Text = '0'
        end
        object RefAntN: TEdit
          Left = 260
          Top = 74
          Width = 51
          Height = 21
          TabOrder = 8
          Text = '0'
        end
        object RefAntU: TEdit
          Left = 312
          Top = 74
          Width = 51
          Height = 21
          TabOrder = 9
          Text = '0'
        end
        object RefPos1: TEdit
          Left = 6
          Top = 36
          Width = 119
          Height = 21
          TabOrder = 1
          Text = '0'
        end
        object RefPos2: TEdit
          Left = 124
          Top = 36
          Width = 119
          Height = 21
          TabOrder = 2
          Text = '0'
        end
        object RefPos3: TEdit
          Left = 244
          Top = 36
          Width = 119
          Height = 21
          TabOrder = 3
          Text = '0'
        end
        object BtnRefPos: TButton
          Left = 345
          Top = 16
          Width = 17
          Height = 17
          Caption = '...'
          Font.Charset = DEFAULT_CHARSET
          Font.Color = clWindowText
          Font.Height = -9
          Font.Name = 'Tahoma'
          Font.Style = []
          ParentFont = False
          TabOrder = 4
          OnClick = BtnRefPosClick
        end
        object RefAntPcv: TCheckBox
          Left = 6
          Top = 58
          Width = 155
          Height = 17
          Caption = 'Antenna Type (*: Auto)'
          TabOrder = 5
          OnClick = RovAntPcvClick
        end
        object RefAnt: TComboBox
          Left = 3
          Top = 74
          Width = 203
          Height = 21
          DropDownCount = 16
          TabOrder = 6
          OnClick = RefAntClick
        end
        object RefPosType: TComboBox
          Left = 6
          Top = 14
          Width = 137
          Height = 21
          Style = csDropDownList
          ItemIndex = 0
          TabOrder = 0
          Text = 'Lat/Lon/Height (deg/m)'
          OnChange = RefPosTypeChange
          Items.Strings = (
            'Lat/Lon/Height (deg/m)'
            'Lat/Lon/Height (dms/m)'
            'X/Y/Z-ECEF (m)'
            'Average of Single Position'
            'Get from Position File'
            'RINEX Header Postion')
        end
        object RefRecTyp: TEdit
          Left = 124
          Top = 98
          Width = 239
          Height = 21
          TabOrder = 10
        end
      end
      object StaPosFile: TEdit
        Left = 3
        Top = 267
        Width = 333
        Height = 21
        TabOrder = 2
      end
      object BtnStaPosFile: TButton
        Left = 351
        Top = 270
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 3
        OnClick = BtnStaPosFileClick
      end
    end
    object TabSheet7: TTabSheet
      Caption = '&Files'
      ImageIndex = 6
      object Label1: TLabel
        Left = 6
        Top = 204
        Width = 102
        Height = 13
        Caption = 'Google Earth Exe File'
      end
      object BtnAntPcvView: TSpeedButton
        Left = 354
        Top = 4
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnAntPcvViewClick
      end
      object Label38: TLabel
        Left = 6
        Top = 2
        Width = 250
        Height = 13
        Caption = 'Satellite/Receiver Antenna PCV File ANTEX/NGS PCV'
      end
      object BtnSatPcvView: TSpeedButton
        Left = 332
        Top = 4
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnSatPcvViewClick
      end
      object Label63: TLabel
        Left = 6
        Top = 60
        Width = 72
        Height = 13
        Caption = 'Geoid Data File'
      end
      object Label15: TLabel
        Left = 6
        Top = 132
        Width = 65
        Height = 13
        Caption = 'DCB Data File'
      end
      object BtnDCBView: TSpeedButton
        Left = 354
        Top = 136
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnDCBViewClick
      end
      object Label47: TLabel
        Left = 6
        Top = 96
        Width = 100
        Height = 13
        Caption = 'Ionosphere Data File'
        Enabled = False
      end
      object BtnIonoView: TSpeedButton
        Left = 354
        Top = 98
        Width = 17
        Height = 15
        Enabled = False
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnIonoViewClick
      end
      object Label51: TLabel
        Left = 6
        Top = 240
        Width = 91
        Height = 13
        Caption = 'BIPM Circular T File'
      end
      object Label5: TLabel
        Left = 6
        Top = 168
        Width = 61
        Height = 13
        Caption = 'ISB Data File'
      end
      object BtnISBView: TSpeedButton
        Left = 354
        Top = 171
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
        OnClick = BtnISBViewClick
      end
      object Label81: TLabel
        Left = 6
        Top = 276
        Width = 65
        Height = 13
        Caption = 'EOP Data File'
      end
      object Label82: TLabel
        Left = 6
        Top = 312
        Width = 60
        Height = 13
        Caption = 'OTL BLQ File'
      end
      object BtnBLQFileView: TSpeedButton
        Left = 355
        Top = 314
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
      end
      object BtnEOPView: TSpeedButton
        Left = 355
        Top = 279
        Width = 17
        Height = 15
        Flat = True
        Glyph.Data = {
          3E020000424D3E0200000000000036000000280000000D0000000D0000000100
          1800000000000802000000000000000000000000000000000000FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080FFFFFFFFFFFFFFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF000000FFFFFF808080808080808080FFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000FFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000FFFFFF00FFFFFF000000
          FFFFFF808080808080808080808080808080808080808080FFFFFF000000FFFF
          FF00FFFFFF000000FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFF000000FFFFFF00FFFFFF00000000000000000000000000000000000000
          0000000000000000000000000000FFFFFF00FFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00FFFFFFFFFFFF
          FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
          FF00}
      end
      object AntPcvFile: TEdit
        Left = 2
        Top = 38
        Width = 353
        Height = 21
        TabOrder = 0
      end
      object BtnAntPcvFile: TButton
        Left = 355
        Top = 40
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 1
        OnClick = BtnAntPcvFileClick
      end
      object BtnGoogleEarthFile: TButton
        Left = 355
        Top = 220
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 5
        OnClick = BtnGoogleEarthFileClick
      end
      object GoogleEarthFile: TEdit
        Left = 2
        Top = 218
        Width = 353
        Height = 21
        TabOrder = 4
      end
      object SatPcvFile: TEdit
        Left = 2
        Top = 16
        Width = 353
        Height = 21
        TabOrder = 6
      end
      object BtnSatPcvFile: TButton
        Left = 355
        Top = 18
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 7
        OnClick = BtnSatPcvFileClick
      end
      object GeoidDataFile: TEdit
        Left = 2
        Top = 74
        Width = 353
        Height = 21
        TabOrder = 2
      end
      object BtnGeoidDataFile: TButton
        Left = 355
        Top = 76
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 3
        OnClick = BtnGeoidDataFileClick
      end
      object DCBFile: TEdit
        Left = 2
        Top = 146
        Width = 353
        Height = 21
        TabOrder = 8
      end
      object BtnDCBFile: TButton
        Left = 355
        Top = 148
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 9
        OnClick = BtnDCBFileClick
      end
      object IonoFile: TEdit
        Left = 2
        Top = 110
        Width = 353
        Height = 21
        Enabled = False
        TabOrder = 10
      end
      object BtnIonoFile: TButton
        Left = 355
        Top = 112
        Width = 17
        Height = 17
        Caption = '...'
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 11
        OnClick = BtnIonoFileClick
      end
      object BIPMCircularTFile: TEdit
        Left = 2
        Top = 254
        Width = 353
        Height = 21
        TabOrder = 12
      end
      object BtnBIPMCircularTFile: TButton
        Left = 355
        Top = 256
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 13
        OnClick = BtnBIPMCircularTFileClick
      end
      object ISBFile: TEdit
        Left = 2
        Top = 182
        Width = 353
        Height = 21
        TabOrder = 14
      end
      object BtnISBFile: TButton
        Left = 355
        Top = 184
        Width = 17
        Height = 17
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 15
        OnClick = BtnISBFileClick
      end
      object EOPFile: TEdit
        Left = 2
        Top = 290
        Width = 353
        Height = 21
        TabOrder = 16
      end
      object BLQFile: TEdit
        Left = 2
        Top = 326
        Width = 353
        Height = 21
        TabOrder = 17
      end
      object BtnBLQFile: TButton
        Left = 355
        Top = 328
        Width = 17
        Height = 19
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 18
      end
      object BtnEOPFile: TButton
        Left = 355
        Top = 292
        Width = 17
        Height = 19
        Caption = '...'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -9
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 19
      end
    end
    object TabSheet6: TTabSheet
      Caption = 'Misc'
      ImageIndex = 7
      object Label19: TLabel
        Left = 170
        Top = 132
        Width = 34
        Height = 13
        Caption = 'Rovers'
      end
      object Label39: TLabel
        Left = -96
        Top = 702
        Width = 37
        Height = 13
        Caption = 'Label39'
      end
      object Label40: TLabel
        Left = 38
        Top = 52
        Width = 150
        Height = 13
        Caption = 'SBAS Satellite Selection (0: All) '
      end
      object Label41: TLabel
        Left = 38
        Top = 8
        Width = 189
        Height = 13
        Caption = 'Time Interpolation of Base Station Data'
      end
      object Label42: TLabel
        Left = 38
        Top = 30
        Width = 121
        Height = 13
        Caption = 'DGPS/DGNSS Corrections'
        Enabled = False
      end
      object Label44: TLabel
        Left = 256
        Top = 132
        Width = 65
        Height = 13
        Caption = 'Base Stations'
      end
      object Label45: TLabel
        Left = 40
        Top = 132
        Width = 67
        Height = 13
        Caption = 'Station ID List'
      end
      object Label85: TLabel
        Left = 36
        Top = 101
        Width = 85
        Height = 13
        Caption = 'RINEX Opt (Base)'
      end
      object Label86: TLabel
        Left = 36
        Top = 79
        Width = 91
        Height = 13
        Caption = 'RINEX Opt (Rover)'
      end
      object RovList: TMemo
        Left = 144
        Top = 146
        Width = 101
        Height = 137
        Lines.Strings = (
          'rover')
        ScrollBars = ssVertical
        TabOrder = 0
      end
      object BaseList: TMemo
        Left = 246
        Top = 146
        Width = 101
        Height = 137
        BevelInner = bvSpace
        Ctl3D = True
        DoubleBuffered = False
        Lines.Strings = (
          'base')
        ParentCtl3D = False
        ParentDoubleBuffered = False
        ParentShowHint = False
        ScrollBars = ssVertical
        ShowHint = False
        TabOrder = 1
      end
      object IntpRefObs: TComboBox
        Left = 232
        Top = 4
        Width = 113
        Height = 21
        Style = csDropDownList
        ItemIndex = 0
        TabOrder = 2
        Text = 'OFF'
        OnChange = FreqChange
        Items.Strings = (
          'OFF'
          'ON')
      end
      object SbasSat: TEdit
        Left = 232
        Top = 48
        Width = 113
        Height = 21
        TabOrder = 3
        Text = '0'
      end
      object ComboBox1: TComboBox
        Left = 232
        Top = 26
        Width = 113
        Height = 21
        Style = csDropDownList
        Enabled = False
        ItemIndex = 0
        TabOrder = 4
        Text = 'SBAS'
        OnChange = FreqChange
        Items.Strings = (
          'SBAS'
          'RTCM')
      end
      object Panel1: TPanel
        Left = 35
        Top = 144
        Width = 103
        Height = 69
        BevelInner = bvRaised
        BevelOuter = bvLowered
        TabOrder = 5
        object SpeedButton1: TSpeedButton
          Left = 4
          Top = 4
          Width = 15
          Height = 17
          Caption = '?'
          Flat = True
          OnClick = SpeedButton1Click
        end
        object Label46: TLabel
          Left = 22
          Top = 4
          Width = 4
          Height = 13
          Caption = ':'
        end
        object Label62: TLabel
          Left = 32
          Top = 5
          Width = 61
          Height = 26
          Caption = 'Keywords in File Path'
          WordWrap = True
        end
        object Label43: TLabel
          Left = 6
          Top = 36
          Width = 20
          Height = 13
          Caption = '#..:'
        end
        object Label65: TLabel
          Left = 32
          Top = 37
          Width = 56
          Height = 26
          Caption = 'Comment in List'
          WordWrap = True
        end
      end
      object RnxOpts1: TEdit
        Left = 144
        Top = 74
        Width = 204
        Height = 21
        TabOrder = 6
      end
      object RnxOpts2: TEdit
        Left = 144
        Top = 96
        Width = 204
        Height = 21
        TabOrder = 7
      end
    end
  end
  object OpenDialog: TOpenDialog
    Filter = 
      'All (*.*)|*.*|PCV File (*.pcv,*.atx)|*.pcv;*.atx|Position File (' +
      '*.pos)|*.pos;*.pos|Options File (*.conf)|*.conf|Execution File (' +
      '*.exe)|*.exe'
    Options = [ofHideReadOnly, ofNoChangeDir, ofEnableSizing]
    OptionsEx = [ofExNoPlacesBar]
    Title = 'Load File'
    Left = 3
    Top = 603
  end
  object SaveDialog: TSaveDialog
    Filter = 'All (*.*)|*.*|Options File (*.conf)|*.conf'
    Options = [ofHideReadOnly, ofNoChangeDir, ofEnableSizing]
    OptionsEx = [ofExNoPlacesBar]
    Title = 'Save File'
    Left = 31
    Top = 603
  end
end
